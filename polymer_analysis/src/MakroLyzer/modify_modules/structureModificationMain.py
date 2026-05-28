from MakroLyzer.input_handling import readXYZ
from MakroLyzer.input_handling import readLMP
from MakroLyzer.input_handling import estimateFrames
from MakroLyzer.input_handling import subgraphSelection
from MakroLyzer import graphs

from MakroLyzer.modify_modules.modifier_registry import MODIFIERS_REGISTRATION

from tqdm import tqdm

def main(args):
    """
    Perfoms the main modification of the polymer structure.
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
        
            
    # Get selections about subgraph selection 
    #################################################################################################
    selection_raw = args.get('subgraphSelection')
    selection_sorted = subgraphSelection.parse_selection_list(selection_raw)
    
    selected_nodes = None
    

    # Build modifiers from registry
    #################################################################################################
    modifiers = {}
    # prepare context for factories
    context = {'boxSize': boxSize}    
    for key, factory in MODIFIERS_REGISTRATION.items():
        try:
            instance = factory(args, **context)
        except Exception:
            instance = None
        if instance is not None:
            modifiers[key] = instance
            if hasattr(instance, 'initialize_output'):
                instance.initialize_output()

    # Main Loop
    ##########################################################################################################################################
    for i, xyz_frame in enumerate(tqdm(read(trajectoryFilePath), total=n_frames, desc="Creating something magical", unit="frame", ncols=100)):
        if i % nthStep != 0:
            continue
         
        boxGraph = graphs.GraphManager(xyz_frame, boxSize=boxSize, vib_factor=vib_factor)
            
        if selection_sorted:
            selected_nodes = boxGraph.select_subgraph_nodes(selection_sorted)
            graph_for_analysis = graphs.GraphManager(boxGraph.subgraph(selected_nodes), boxSize=boxSize, vib_factor=vib_factor)
            boxGraph = graph_for_analysis
        
        # Run all registered analyzers for this frame with the graph_for_analyis
        for modifier in modifiers.values():
            modifier.run(boxGraph, i)
    ##########################################################################################################################################
