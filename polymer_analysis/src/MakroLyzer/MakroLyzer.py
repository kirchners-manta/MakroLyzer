import sys

from MakroLyzer.input_handling import inputHandlingMain
from MakroLyzer.structure_modules import structureAnalysisMain
from MakroLyzer.modify_modules import structureModificationMain

def main():
    """
    Main function to run the PolyLyzer program.
    """
    # Get command line arguments and xyz data
    args = inputHandlingMain.main(sys.argv)
    
    # Call the structure analysis of the polymer structure
    structureAnalysisMain.main(args)
    
    # Call the modify modules 
    structureModificationMain.main(args)
    
if __name__ == "__main__":
    main()