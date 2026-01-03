from MakroLyzer.input_handling import readXYZ
from MakroLyzer.input_handling import readLMP
from MakroLyzer.input_handling import estimateFrames
from MakroLyzer import graphs
from MakroLyzer.structure_modules import readPatterns

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
        'dihedral_list': [],
        'cisTrans': [],
        'AARamachandran': [],
        'Rg': [],
        'hbonds': [],
        'subgraph_coords': [],
        'anisotropy_factor': [],
        'asphericity_parameter': [],
        'orderParameter': [],
        'RingStrandCount': [],

        # Output file names
        'formulas_file': args['formula_file'],
        'noSub_file': args['noSub_file'],
        'distances_file': args['e2e_file'],
        'dihedrals_file': args['dihedral_file'],
        'cisTrans_file': args['CisTrans_file'],
        'AARamachandran_file': args['AARamachandran_file'],
        'Rg_file': args['Rg_file'],
        'hbonds_file': args['hbonds_file'],
        'subgraph_coords_file': args['subgraph_coord_file'],
        'anisotropy_file': args['anisotropy_file'],
        'asphericity_file': args['asphericity_file'],
        'orderParameter_file': args['order_file'],
        'RingStrandCount_file': args['RingStrandCount_file']
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
            
            
        # Subgraph coordinates
        if args['subgraph_coords']:
            results['subgraph_coords'].append(get_subgraph_coords(boxGraph))
    return results