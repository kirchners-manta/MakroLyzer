from MakroLyzer.input_handling import readXYZ
from MakroLyzer.input_handling import readLMP
from MakroLyzer.input_handling import estimateFrames
from MakroLyzer.input_handling import subgraphSelection
from MakroLyzer import graphs

from MakroLyzer.structure_modules.analyzer_registry import ANALYZERS_REGISTRATION
from MakroLyzer.structure_modules.backbone_cache import BackboneCache

from tqdm import tqdm


def main(args):
    """
    Perfoms the main analysis of the polymer structure.
    """
    
    # Get the trajectory file path
    #################################################################################################
    if args['xyzFile']:
        trajectoryFilePath = args['xyzFile']
        n_frames = estimateFrames.EstimateFrames.estimateFramesXYZ(trajectoryFilePath)
        read = readXYZ.readXYZ
    elif args['lmpFile']:
        trajectoryFilePath = args['lmpFile']
        n_frames = estimateFrames.EstimateFrames.estimateFramesLMP(trajectoryFilePath)
        read = readLMP.readLMP
        
        
    # Get the modulo for reading frames 
    #################################################################################################
    nthStep = args.get('nthStep', 1)
    
    
    # Get the box size 
    #################################################################################################
    boxSize = args.get('BoxSize', None)
    vib_factor = args.get('vibFactor', None)
        # check if the order parameter is provided
    if args['orderParameter']:
        BoxSize, NoCellsPerDim, MolecularVectorLength = args['orderParameter']
        op_vector_source = args.get('op_vector_source', 'atom')
        if op_vector_source == 'ring':
            op_vector_source = 'ring-center'
        op_ring_cycle_size = args.get('op_ring_size')
        if op_vector_source == 'ring-normal':
            MolecularVectorLength = None
    if not boxSize and args['orderParameter']:
        boxSize = BoxSize[0]
        
            
    # Get selections about static-topology and subgraph selection 
    #################################################################################################
    static_topology =args.get('staticTopology', False)
    selection_raw = args.get('subgraphSelection')
    selection_sorted = subgraphSelection.parse_selection_list(selection_raw)
    
    base_graph = None
    selected_nodes = None
    
    
    # Build analyzers from registry
    #################################################################################################
    analyzers = {}
    analyzer_graph_routes = []
    # prepare context for factories
    context = {'boxSize': boxSize, 'static_topology': static_topology}
    if args.get('orderParameter'):
        context.update(
            {
                'BoxSize': BoxSize,
                'NoCellsPerDim': NoCellsPerDim,
                'MolecularVectorLength': MolecularVectorLength,
                'op_vector_source': op_vector_source,
                'op_ring_size': op_ring_cycle_size,
            }
        )
    if static_topology and (
        args.get('orderParameter') or args.get('dihedral') or args.get('endToEndDistance')
    ):
        # We only need to cache the backbone if the user wants a static-topology analysis and at least one of the analyzers 
        # that uses the backbone is selected
        context['backbone_cache'] = BackboneCache()

    for key, factory in ANALYZERS_REGISTRATION.items():
        try:
            instance = factory(args, **context)
        except Exception as exc:
            print(f"Analyzer '{key}' failed to initialize: {exc}")
            raise
        if instance is not None:
            analyzers[key] = instance
            analyzer_graph_routes.append((instance, getattr(instance, 'requires_full_graph', False)))
            if hasattr(instance, 'initialize_output'):
                instance.initialize_output()


    # Main Loop
    ##########################################################################################################################################
    for i, xyz_frame in enumerate(tqdm(read(trajectoryFilePath), total=n_frames, desc="Creating something magical", unit="frame", ncols=100)):
        if i % nthStep != 0:
            continue
         
        # The base_graph contains all atoms as nodes that we read from the specific file, 
        # not only the user selection.
        if base_graph is None or not static_topology:
            # We dont have a base_graph OR the user does not want a static-topology analyis
            # The if statement is only False if base_graph is not None and static_topology is True 
            # (skipped when: We already have a base graph and we want static_topology)
            # -> We need to (re)build the graph
            boxGraph = graphs.GraphManager(xyz_frame, boxSize=boxSize, vib_factor=vib_factor)
            
            if static_topology:
                # The user wants a static-topology analyis
                # We only arive here if we didn't have a base_graph before
                # -> We save the generated graph as base_graph
                base_graph = boxGraph
                
                if selection_sorted:
                    # The used additionally selected specific subgraphs
                    # We only arive here if we didn't have a base_graph before and the user wants static_topology
                    # -> We need to select the nodes in the graph
                    selected_nodes = boxGraph.select_subgraph_nodes(selection_sorted)
                    
        else:
            # We already have a base graph and we want static_topology
            # -> we reuse the base_graph from the iteration before and update the coordinates
            boxGraph = base_graph
            boxGraph.update_coordinates(xyz_frame, boxSize=boxSize)
            
        # -> Now we select the graph to perform the analysis with
        graph_for_analysis = boxGraph
        
        if selection_sorted:
            # If the user selected specific subgraphs, we need to adjust the graph_for_analysis
            
            if selected_nodes is None or not static_topology:
                # The nodes are not selected yet, or they are selected, but the user does not want a static topology
                # (skipped when: The nodes are selected and the user wants static topology analysis)
                selected_nodes = boxGraph.select_subgraph_nodes(selection_sorted)

            # -> We need to update the graph_for_analyis 
            # We only analyze the selected subgraphs (They are all contained in boxGraph.subgraph(selected_nodes))
            graph_for_analysis = graphs.GraphManager(boxGraph.subgraph(selected_nodes), boxSize=boxSize, vib_factor=vib_factor)
        
        # Run all registered analyzers for this frame. Full-graph analyses must
        # not change the graph passed to analyzers that respect global -sel.
        for analyzer, requires_full_graph in analyzer_graph_routes:
            if requires_full_graph:
                analyzer_graph = boxGraph
            else:
                analyzer_graph = graph_for_analysis
            analyzer.run(analyzer_graph, i)
    ##########################################################################################################################################
    
    
    # Finalize output ('selcetion'-mode) 
    #################################################################################################
    for analyzer_name, analyzer in analyzers.items():
        analyzer.finalize_output()
