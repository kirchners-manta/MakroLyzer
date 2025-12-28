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

from MakroLyzer.structure_modules.endToEndDistance import end_to_end_dist
from MakroLyzer.structure_modules.dihedrals import get_all_dihedrals, get_CisTrans, get_Ramachandran
from MakroLyzer.structure_modules.subgraphCoords import get_subgraph_coords
from MakroLyzer.structure_modules.orderParameter import get_order_parameter, get_S_from_Q

from tqdm import tqdm

def main(args):
    """
    Perfoms the main analysis of the polymer structure.
    """
    
    
    # Vorrübergende Lösung
    analyzers = {}
    
    if args['hydrogenBonds']:
        cutoffs = args['hydrogenBonds']
        output_handler = OutputHandler(args['hbonds_file'], mode='streaming')
        analyzers['hbonds'] = HBondsAnalyzer(cutoffs, output_handler)
        analyzers['hbonds'].initialize_output()
        
    if args['anisotropyFactor']:
        output_handler = OutputHandler(args['anisotropy_file'], mode='collect')
        analyzers['anisotropy'] = AnisotropyAnalyzer(output_handler)
        analyzers['anisotropy'].initialize_output()
        
    if args['asphericityParameter']:
        output_handler = OutputHandler(args['asphericity_file'], mode='streaming')
        analyzers['asphericity'] = AsphericityAnalyzer(output_handler)
        analyzers['asphericity'].initialize_output()
        
    if args['radiusOfGyration']:
        output_handler = OutputHandler(args['Rg_file'], mode='collect')
        analyzers['radiusOfGyration'] = RadiusOfGyrationAnalyzer(output_handler)
        analyzers['radiusOfGyration'].initialize_output()
        
    if args['noSubgraphs']:
        output_handler = OutputHandler(args['noSub_file'], mode='streaming')
        analyzers['noSubgraphs'] = MoleculeCountAnalyzer(output_handler)
        analyzers['noSubgraphs'].initialize_output()
    
    # create empty lists to store results
    results = {
        'formulas': [],
        'distances': [],
        'dihedrals': [],
        'dihedral_list': [],
        'cisTrans': [],
        'AARamachandran': [],
        'subgraph_coords': [],
        'orderParameter': [],

        # Output file names
        'formulas_file': args['formula_file'],
        'distances_file': args['e2e_file'],
        'dihedrals_file': args['dihedral_file'],
        'cisTrans_file': args['CisTrans_file'],
        'AARamachandran_file': args['AARamachandran_file'],
        'subgraph_coords_file': args['subgraph_coord_file'],
        'orderParameter_file': args['order_file']
    }
    
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
            
        # End-to-end distance     
        if args['endToEndDistance']:
            dist = end_to_end_dist(boxGraph)
            results.setdefault('distances', []).append(dist)
          
        # Dihedral angles      
        if args['dihedral']:
            if args['dihedral_range'] == 'abs':
                results['dihedrals'].append(get_all_dihedrals(boxGraph, sign=None)[0])
                results['dihedral_list'].append(get_all_dihedrals(boxGraph, sign=None)[1])
            elif args['dihedral_range'] == 'nonabs':
                results['dihedrals'].append(get_all_dihedrals(boxGraph, sign=True)[0])
                results['dihedral_list'].append(get_all_dihedrals(boxGraph, sign=True)[1])

        # Cis Trans counts
        if args['cisTrans']:
            results['cisTrans'].append(get_CisTrans(boxGraph))

        # Amino Acid Ramachandran plot data
        if args['AminoAcidRamachandran']:
            results['AARamachandran'].append(get_Ramachandran(boxGraph))
            
        # Subgraph coordinates
        if args['subgraph_coords']:
            results['subgraph_coords'].append(get_subgraph_coords(boxGraph))            
            
        # Order parameter
        if args['orderParameter']:
            boxSize, n, unitSize = args['orderParameter']
            results['orderParameter'].append(get_S_from_Q(
                boxGraph, boxSize, n, unitSize))
            
    # Finalize output for class-based analyzers ---------------------------------
    for analyzer_name, analyzer in analyzers.items():
        analyzer.finalize_output()

    return results