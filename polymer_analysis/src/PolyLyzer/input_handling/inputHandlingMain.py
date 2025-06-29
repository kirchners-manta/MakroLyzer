import sys

from PolyLyzer.input_handling import readInput
from PolyLyzer.input_handling import checkInput

def main(args):
    """
    Input reading and checking.
    """    
    
    args = readInput.readCommandLine()
    
    try:
        checkInput.checkInput(args)
    except checkInput.FileNotFoundError as e:
        print(e)
        sys.exit(1)
    except checkInput.InvalidFileFormatError as e:
        print(e)
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)
        
    return args