from PolyLyzer.input_handling import readXYZ
from PolyLyzer.input_handling import readLMP
from PolyLyzer.input_handling import estimateFrames
from PolyLyzer.structure_modules import graphs
from PolyLyzer.structure_modules import readPatterns
from PolyLyzer.structure_modules.endToEndDistance import end_to_end_dist
from PolyLyzer.structure_modules.countSubgraphs import count_subgraphs
from PolyLyzer.structure_modules.dihedrals import get_all_dihedrals, get_CisTrans
from PolyLyzer.structure_modules.radiusOfGyration import get_radius_of_gyration
from PolyLyzer.structure_modules.anisotropy import get_anisotropy_factor
from PolyLyzer.structure_modules.asphericityParameter import get_asphericity_parameter
from PolyLyzer.structure_modules.hbonds import get_Hbonds
from PolyLyzer.structure_modules.subgraphCoords import get_subgraph_coords
from PolyLyzer.structure_modules.orderParameter import get_order_parameter

from tqdm import tqdm

def main(args):
    """
    Perfoms the main analysis of the polymer structure.
    """
    
    # create empty lists to store results
    results = {
        'formulas': [],
        'noSub': [],
        'distances': [],
        'dihedrals': [],
        'cisTrans': [],
        'Rg': [],
        'hbonds': [],
        'subgraph_coords': [],
        'anisotropy_factor': [],
        'asphericity_parameter': [],
        'orderParameter': [],
        
        # Output file names
        'formulas_file': args['formula_file'],
        'noSub_file': args['noSub_file'],
        'distances_file': args['e2e_file'],
        'dihedrals_file': args['dihedral_file'],
        'cisTrans_file': args['CisTrans_file'],
        'Rg_file': args['Rg_file'],
        'hbonds_file': args['hbonds_file'],
        'subgraph_coords_file': args['subgraph_coord_file'],
        'anisotropy_file': args['anisotropy_file'],
        'asphericity_file': args['asphericity_file'],
        'orderParameter_file': args['order_file']
    }
    
    # Get the trajectory file path 
    if args['xyzFile']:
        trajectoryFilePath = args['xyzFile']
        n_frames = estimateFrames.EstimateFrames.estimateFramesXYZ(trajectoryFilePath)
        read = readXYZ.readXYZ
    elif args['lmpFile']:
        trajectoryFilePath = args['lmpFile']
        n_frames = estimateFrames.EstimateFrames.estimateFrameLMP(trajectoryFilePath)
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
            
    
    for i, xyz_frame in enumerate(tqdm(read(trajectoryFilePath),total=n_frames, desc="Creating something magical", unit="frame", ncols=100)):
        if i % nthStep != 0:
            continue
        
        # Get Graph object of the polymer box
        boxGraph = graphs.GraphManager(xyz_frame, boxSize=boxSize)
        
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
            
        # Number of subgraphs
        if args['noSubgraphs']:
            noSub = count_subgraphs(boxGraph)
            results['noSub'].append(noSub)
           
        # End-to-end distance     
        if args['endToEndDistance']:
            dist = end_to_end_dist(boxGraph)
            results.setdefault('distances', []).append(dist)
          
        # Dihedral angles      
        if args['dihedral']:
            if args['dihedral_range'] == 'abs':
                results['dihedrals'].append(get_all_dihedrals(boxGraph, sign=None))
            elif args['dihedral_range'] == 'nonabs':
                results['dihedrals'].append(get_all_dihedrals(boxGraph, sign=True))
                
        # Cis Trans counts
        if args['cisTrans']:
            results['cisTrans'].append(get_CisTrans(boxGraph))
           
        # Radius of gyration 
        if args['radiusOfGyration']:
            results['Rg'].append(get_radius_of_gyration(boxGraph))
          
        # Hydrogen bonds  
        if args['hydrogenBonds']:
            cutoffs = args['hydrogenBonds']
            results['hbonds'].append(get_Hbonds(boxGraph, cutoffs))
            
        # Subgraph coordinates
        if args['subgraph_coords']:
            results['subgraph_coords'].append(get_subgraph_coords(boxGraph))
            
        # Anisotropy factor
        if args['anisotropyFactor']:
            results['anisotropy_factor'].append(get_anisotropy_factor(boxGraph))
            
            
        # Asphericity parameter
        if args['asphericityParameter']:
            results['asphericity_parameter'].append(get_asphericity_parameter(boxGraph))
            
        # Order parameter
        if args['orderParameter']:
            boxSize, n, unitSize = args['orderParameter']
            results['orderParameter'].append(get_order_parameter(
                boxGraph, boxSize, n, unitSize))

    return results