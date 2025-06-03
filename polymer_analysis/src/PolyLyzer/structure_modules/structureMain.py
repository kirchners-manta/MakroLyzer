import sys
from IPython.display import display

from PolyLyzer.structure_modules import graphs
from PolyLyzer.structure_modules import readPatterns
from PolyLyzer.structure_modules.endToEndDistance import end_to_end_dist
from PolyLyzer.structure_modules.dihedrals import get_all_dihedrals, get_CisTrans

def main(args, xyz):
    """
    Perfoms the main analysis of the polymer structure.
    """
    
    # Get Graph object of the polymer box
    boxGraph = graphs.GraphManager(xyz)
    
    # Repeating units for the Polymer
    if args['patternFile']:
        patternInfo = readPatterns.readPattern(args['patternFile'])
        boxGraph.find_and_tag_patterns(patternInfo['pattern'].values[0], patternInfo['element'].values[0])
        boxGraph.write_fragment_data_to_csv("repeating_units.csv")
        # display boxGraph
        #for node in boxGraph.nodes(data=True):
        #    display(node)
        
    # Saturation
    if args['saturation']:
        boxGraph.saturate()
        boxGraph.write_xyz("n66_saturated.xyz")
        formulas = boxGraph.get_chemicalFormulas()
        print("Saturation formulas:")
        for formula in formulas:
            print(formula)
            
    # End-to-end distance
    if args['endToEndDistance']:
        distances = end_to_end_dist(boxGraph)
        print("End-to-end distances:")
        for distance in distances:
            print(distance)
            
    # Dihedral angles
    if args['dihedral']:
        # test if --dihedral-range and --dihedral-file are provided
        if args['dihedral_range'] == 'abs':
            dihedrals = get_all_dihedrals(boxGraph, file=args['dihedral_file'], sign=None)
        elif args['dihedral_range'] == 'nonabs':
            dihedrals = get_all_dihedrals(boxGraph,  file=args['dihedral_file'], sign=True)
            
    # Cis and Trans counts
    if args['cisTrans']:
        ct = get_CisTrans(boxGraph, file=args['CisTrans_file'])