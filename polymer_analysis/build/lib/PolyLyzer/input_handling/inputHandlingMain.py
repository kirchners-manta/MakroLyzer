import sys

from PolyLyzer.input_handling import readInput
from PolyLyzer.input_handling import checkInput
from IPython.display import display

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
        
    xyzFilePath = args['xyzFile']
    xyz = readInput.readXYZ(xyzFilePath)
    
    # Wrap the XYZ coordinates if periodic boundary conditions are applied
    if args['PBC_xyz']:
        x = args['PBC_xyz'][0]
        y = args['PBC_xyz'][1]
        z = args['PBC_xyz'][2]
        xyz = readInput.wrapXYZ(xyz, x, y, z)
        
    return args, xyz