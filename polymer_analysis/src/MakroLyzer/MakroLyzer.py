import sys

from MakroLyzer.input_handling import inputHandlingMain
from MakroLyzer.structure_modules import structureAnalysisMain
from MakroLyzer.modify_modules import structureModificationMain

def main():
    """
    Main function to run the MakroLyzer program.
    """
    # Get command line arguments and xyz data
    analyzer, modifier,args = inputHandlingMain.main(sys.argv)
        
    # Call the structure analysis of the polymer structure
    if analyzer:
        structureAnalysisMain.main(args)
    
    # Call the modify modules 
    if modifier:
        structureModificationMain.main(args)
    
if __name__ == "__main__":
    main()