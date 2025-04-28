import sys
from IPython.display import display

from PolyLyzer.input_handling import inputHandlingMain
from PolyLyzer.structure_modules import graphs
from PolyLyzer.structure_modules import structureMain


def main():
    """
    Main function to run the PolyLyzer program.
    """
    
    # Get command line arguments and xyz data
    args, xyz = inputHandlingMain.main(sys.argv)
    
    # Get Graph object of the polymer box
    if args['PBC_xyz']:
        print("Periodic boundary conditions applied: ", args['PBC_xyz'])
        boxGraph = graphs.GraphManager(xyz, args['PBC_xyz'])
    else:
        boxGraph = graphs.GraphManager(xyz)
    
    # Call the main analysis of the polymer structure
    structureMain.main(args, boxGraph)
    
    # End to end distance
    #print("End to end distance: ", boxGraph.end_to_end_dist().mean(), " +/- ", boxGraph.end_to_end_dist().std())
    #print("Ensemble average end to end distance: ", boxGraph.end_to_end_dist_ensemble())

    

    
if __name__ == "__main__":
    main()