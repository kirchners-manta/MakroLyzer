import sys

from MakroLyzer.input_handling import readInput
from MakroLyzer.input_handling import checkInput

def main(args):
    """
    Input reading and checking.
    """    
    
    args = readInput.readCommandLine()
    analyzer = True
    modifier = True
    
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
        
    # Check if args from analyzer_group are called
    analyzer_args = [
        'hydrogenBonds', 'anisotropyFactor', 'asphericityParameter',
        'radiusOfGyration', 'NoSubgraphs', 'endToEndDistance', 'orderParameter',
        'formula', 'Ramachandran', 'dihedral'
    ]
    if not any(arg in args and args[arg] is not None for arg in analyzer_args):
        analyzer = False
        
    # Check if args from modifier_group are called
    modifier_args = [
        'functionalizePE', 'functionalizePEsurface', 'patternFile', 'saturation', 'subgraph_coords',
        'wrapGraphPrint'
    ]
    if not any(arg in args and args[arg] is not None for arg in modifier_args):
        modifier = False

    return analyzer, modifier, args
