from MakroLyzer.input_handling import readXYZ
from MakroLyzer.input_handling import readLMP
from MakroLyzer.input_handling import estimateFrames
from MakroLyzer.structure_modules import graphs
from MakroLyzer.structure_modules import readPatterns
from MakroLyzer.structure_modules.structureBase import OutputHandler

from MakroLyzer.structure_modules.Hbonds import HBondsAnalyzer
from MakroLyzer.structure_modules.Anisotropy import AnisotropyAnalyzer
from MakroLyzer.structure_modules.Asphericity import AsphericityAnalyzer
from MakroLyzer.structure_modules.RadiusOfGyration import RadiusOfGyrationAnalyzer
from MakroLyzer.structure_modules.MoleculeCount import MoleculeCountAnalyzer
from MakroLyzer.structure_modules.EndToEndDistance import EndToEndDistanceAnalyzer
from MakroLyzer.structure_modules.OrderParameter import OrderParameterAnalyzer
from MakroLyzer.structure_modules.Ramachandran import RamachandranAnalyzer
from MakroLyzer.structure_modules.Dihedrals import DihedralsAnalyzer

from MakroLyzer.structure_modules.subgraphCoords import get_subgraph_coords

from tqdm import tqdm

def main(args):
    """
    Perfoms the main analysis of the polymer structure.
    """
    # Get the trajectory file path 
    if args['xyzFile']:
        trajectoryFilePath = args['xyzFile']
        n_frames = estimateFrames.EstimateFrames.estimateFramesXYZ(trajectoryFilePath)
        read = readXYZ.readXYZ
    elif args['lmpFile']:
        trajectoryFilePath = args['lmpFile']
        n_frames = estimateFrames.EstimateFrames.estimateFramesLMP(trajectoryFilePath)
        read = readLMP.readLMP
        
    # Get the modulo for reading frames
    nthStep = args.get('nthStep', 1)
    
    # Get the box size
    boxSize = args.get('BoxSize', None)
    if not boxSize:
        # check if the order parameter is provided
        if args['orderParameter']:
            BoxSize, n, unitSize = args['orderParameter']
            boxSize = BoxSize[0]
            
    
    # Vorrübergende Lösung -----------------------------------------------------------------------
    analyzers = {}
    
    if args['hydrogenBonds']:
        cutoffs = args['hydrogenBonds']
        output_handler = OutputHandler(args['hbonds_file'], mode='streaming')
        analyzers['hbonds'] = HBondsAnalyzer(cutoffs, output_handler)
        analyzers['hbonds'].initialize_output()
        
    if args['anisotropyFactor']:
        output_handler = OutputHandler(args['anisotropy_file'], mode='streaming')
        analyzers['anisotropy'] = AnisotropyAnalyzer(output_handler)
        analyzers['anisotropy'].initialize_output()
        
    if args['asphericityParameter']:
        output_handler = OutputHandler(args['asphericity_file'], mode='streaming')
        analyzers['asphericity'] = AsphericityAnalyzer(output_handler)
        analyzers['asphericity'].initialize_output()
        
    if args['radiusOfGyration']:
        output_handler = OutputHandler(args['Rg_file'], mode='streaming')
        analyzers['radiusOfGyration'] = RadiusOfGyrationAnalyzer(output_handler)
        analyzers['radiusOfGyration'].initialize_output()
        
    if args['noSubgraphs']:
        output_handler = OutputHandler(args['noSub_file'], mode='streaming')
        analyzers['noSubgraphs'] = MoleculeCountAnalyzer(output_handler)
        analyzers['noSubgraphs'].initialize_output()
    
    if args['endToEndDistance']:
        output_handler = OutputHandler(args['e2e_file'], mode='streaming')
        analyzers['endToEndDistance'] = EndToEndDistanceAnalyzer(output_handler)
        analyzers['endToEndDistance'].initialize_output()
        
    if args['orderParameter']:
        output_handler = OutputHandler(args['order_file'], mode='streaming')
        analyzers['orderParameter'] = OrderParameterAnalyzer(boxSize, n, unitSize, output_handler)
        analyzers['orderParameter'].initialize_output()
        
    if args['AminoAcidRamachandran']:
        output_handler = OutputHandler(args['AARamachandran_file'], mode='streaming')
        analyzers['AminoAcidRamachandran'] = RamachandranAnalyzer(output_handler)
        analyzers['AminoAcidRamachandran'].initialize_output()
        
    if args['dihedral'] or args['cisTrans']:
        dihedral_handler = OutputHandler(args['dihedral_file'], mode='collect') if args['dihedral'] else None
        dihedral_list_handler = OutputHandler(args['dihedral_file'].replace('.csv', '_list.csv'), mode='streaming') if args['dihedral'] else None
        cistrans_handler = OutputHandler(args['CisTrans_file'], mode='streaming') if args['dihedral'] else None
        analyzers['dihedrals'] = DihedralsAnalyzer(
            dihedral_output_handler=dihedral_handler,
            dihedral_list_output_handler=dihedral_list_handler,
            cistrans_output_handler=cistrans_handler,
            dihedral_range=args['dihedral_range']
        )
        analyzers['dihedrals'].initialize_output()
        
    # ------------------------------------------------------------------------------------------------
    
    # create empty lists to store results
    results = {
        'formulas': [],
        'subgraph_coords': [],

        # Output file names
        'formulas_file': args['formula_file'],
        'subgraph_coords_file': args['subgraph_coord_file']
    }
    

    # Loop ----------------------------------------------------
    for i, xyz_frame in enumerate(tqdm(read(trajectoryFilePath),total=n_frames, desc="Creating something magical", unit="frame", ncols=100)):
        if i % nthStep != 0:
            continue
        
        # Get Graph object of the polymer box
        boxGraph = graphs.GraphManager(xyz_frame, boxSize=boxSize)
        
        # ---------------------------------------------------------------------
        # Hydrogen bonds  
        if 'hbonds' in analyzers:
            analyzers['hbonds'].run(boxGraph, i)
            
        # Anisotropy factor
        if 'anisotropy' in analyzers:
            analyzers['anisotropy'].run(boxGraph, i)  
            
        # Asphericity parameter
        if 'asphericity' in analyzers:
            analyzers['asphericity'].run(boxGraph, i) 
            
        # Radius of Gyration
        if 'radiusOfGyration' in analyzers:
            analyzers['radiusOfGyration'].run(boxGraph, i)
            
        # Molecule Count
        if 'noSubgraphs' in analyzers:
            analyzers['noSubgraphs'].run(boxGraph, i)
            
        # EndToEnd Distances
        if 'endToEndDistance' in analyzers:
            analyzers['endToEndDistance'].run(boxGraph, i)
            
        # Order Parameter
        if 'orderParameter' in analyzers:
            analyzers['orderParameter'].run(boxGraph, i)
            
        # Ramachandran
        if 'AminoAcidRamachandran' in analyzers:
            analyzers['AminoAcidRamachandran'].run(boxGraph, i)
            
        # Dihedrals and Cis-Trans
        if 'dihedrals' in analyzers:
            analyzers['dihedrals'].run(boxGraph, i)
            
        # ---------------------------------------------------------------------
        
        # Repeating units for the Polymer
        if args['patternFile']:
            patternInfo = readPatterns.readPattern(args['patternFile'])
            boxGraph.find_and_tag_patterns(
                patternInfo['pattern'].values[0],
                patternInfo['element'].values[0]
            )
            # Include frame number in the output file
            repeating_units_file = f"{args['repeatingUnits_file'].rsplit('.', 1)[0]}_frame_{i}.csv"
            boxGraph.write_fragment_data_to_csv(repeating_units_file)

            
        # Saturation
        if args['saturation']:
            boxGraph.saturate()
            saturation_file = f"{args['saturation_file'].rsplit('.', 1)[0]}_frame_{i}.xyz"
            boxGraph.write_xyz(saturation_file)
            
        # Chemical formulas
        if args['formula']:
            formulas = boxGraph.get_chemicalFormulas()
            results['formulas'].append(formulas)
            
        # Subgraph coordinates
        if args['subgraph_coords']:
            results['subgraph_coords'].append(get_subgraph_coords(boxGraph))            
            
    # Finalize output for class-based analyzers ---------------------------------
    for analyzer_name, analyzer in analyzers.items():
        analyzer.finalize_output()

    return results