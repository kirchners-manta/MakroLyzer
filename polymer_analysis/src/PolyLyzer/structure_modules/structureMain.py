from PolyLyzer.input_handling import readXYZ
from PolyLyzer.structure_modules import graphs
from PolyLyzer.structure_modules import readPatterns
from PolyLyzer.structure_modules.endToEndDistance import end_to_end_dist
from PolyLyzer.structure_modules.dihedrals import get_all_dihedrals, get_CisTrans
from PolyLyzer.structure_modules.radiusOfGyration import get_radius_of_gyration
from PolyLyzer.structure_modules.hbonds import get_Hbonds

def main(args):
    """
    Perfoms the main analysis of the polymer structure.
    """
    
    # Get the trajectory file path
    trajectoryFilePath = args['xyzFile']
    
    # create empty lists to store results
    results = {
        'formulas': [],
        'distances': [],
        'dihedrals': [],
        'cisTrans': [],
        'Rg': [],
        'hbonds': [],
        
        # Output file names
        'formulas_file': args['formula_file'],
        'distances_file': args['e2e_file'],
        'dihedrals_file': args['dihedral_file'],
        'cisTrans_file': args['CisTrans_file'],
        'Rg_file': args['Rg_file'],
        'hbonds_file': args['hbonds_file']
    }
    
    
    for i, xyz_frame in enumerate(readXYZ.readXYZ(trajectoryFilePath)):
        print(f"Processing frame {i}")
        # Get Graph object of the polymer box
        boxGraph = graphs.GraphManager(xyz_frame)
        
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
            boxGraph.write_xyz(args['saturation_file'])
            
        # Chemical formulas
        if args['formula']:
            formulas = boxGraph.get_chemicalFormulas()
            results['formulas'].append(formulas)
           
        # End-to-end distance     
        if args['endToEndDistance']:
            results['distances'].append(end_to_end_dist(boxGraph))
          
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
            TypesDistances = args['hydrogenBonds']
            results['hbonds'].append(get_Hbonds(boxGraph, TypesDistances))
            

    return results