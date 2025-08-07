import sys
from IPython.display import display

from MakroLyzer.input_handling import inputHandlingMain
from MakroLyzer.structure_modules import graphs
from MakroLyzer.structure_modules import structureMain
from MakroLyzer.output_handling import outputHandlingMain


def main():
    """
    Main function to run the PolyLyzer program.
    """
    # Get command line arguments and xyz data
    args = inputHandlingMain.main(sys.argv)
    
    # Call the main analysis of the polymer structure
    results = structureMain.main(args)
    # Call the output handling module to save results
    outputHandlingMain.main(results)

    
if __name__ == "__main__":
    main()