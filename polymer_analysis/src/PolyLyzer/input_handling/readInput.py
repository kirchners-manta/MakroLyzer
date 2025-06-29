import argparse
import sys
import pandas as pd

def readCommandLine() -> dict:
    
    """
    Read input parameters from command line arguments.
    
    Returns:
        dict: Dictionary containing the input parameters.
    """
    
    # Create the argument parser
    parser = argparse.ArgumentParser(description='PolyLyzer - Polymer Analysis Tool')
    
    # Add arguments
    parser.add_argument('-xyz', '--xyzFile', 
                        help='Path to the XYZ/trajectory file', 
                        required=True) 
    
    parser.add_argument('-p', '--patternFile', 
                        help='Path to the TXT file -> Finds repeating units (default: false)'
    )
    
    parser.add_argument('--repeatingUnits-file',
                        help='Output file name for repeating units (default: repeatingUnits.csv)',
                        default='repeatingUnits.csv'
    )
    
    parser.add_argument('-s', '--saturation', 
                        help='Saturate the ends polymers (default: false)', 
                        action='store_true'
    )
    
    parser.add_argument('--saturation-file',
                        help='Output file name for saturated polymers (default: saturatedPolymers.xyz)',
                        default='saturatedPolymers.xyz'
    )
    
    parser.add_argument('-f', '--formula',
                        help='Get chemical formulas of the polymer (default: false)', 
                        action='store_true'
    )
    parser.add_argument('--formula-file',
                        help='Output file name for chemical formulas (default: chemicalFormulas.csv)',
                        default='chemicalFormulas.csv'
    )
    
    parser.add_argument('-e2e', '--endToEndDistance', 
                        help='Calculate end-to-end distance (default: false)', 
                        action='store_true'
    )
    
    parser.add_argument(
                        '--e2e-file',
                        help='Output file name for end-to-end distances (default: endToEndDistances.csv)',
                        default='endToEndDistances.csv'
    )
    
    parser.add_argument(
                        '-d', '--dihedral',
                        help='Calculate dihedral angles (default: false)',
                        action='store_true'
    )
    parser.add_argument(
                        '--dihedral-range',
                        help='Range of dihedral angles - abs (0-180), nonabs (-180-180) (default: 0-180)',
                        choices=['abs', 'nonabs'],
                        default='abs'
    )
    parser.add_argument(
                        '--dihedral-file',
                        help='Output file name (default: dihedrals.csv)',
                        default='dihedrals.csv'
    )
    parser.add_argument(
                        '-ct', '--cisTrans',
                        help='Calculate cis and trans counts (default: false)',
                        action='store_true'
    )
    parser.add_argument(
                        '--CisTrans-file',
                        help='Output file name (default: CisTrans.csv)',
                        default='CisTrans.csv'
    )
    parser.add_argument(
                        '-r', '--radiusOfGyration',
                        help='Calculate radius of gyration (default: false)',
                        action='store_true'
    )
    parser.add_argument(
                        '--Rg-file',
                        help='Output file name (default: radiusOfGyration.csv)',
                        default='radiusOfGyration.csv'
    )
    parser.add_argument(
                        '-hb', '--hydrogenBonds',
                        nargs='+',
                        help="List of (element:distance) tuples for hydrogen bonds (e.g., -hb O:3.5 N:2.8)",
                        type=element_distance_tuple
    )
    
    parser.add_argument(
                        '--hbonds-file',
                        help='Output file name for hydrogen bonds (default: hydrogenBonds.csv)',
                        default='hydrogenBonds.csv'
    )
    
    parser.add_argument(
                        '-sub', '--subgraph-coords',
                        help='Get subgraph-coordinates (default: false)',
                        action='store_true'   
    )
    
    parser.add_argument(
                        '--subgraph-coord-file',
                        help='Output file name for subgraph coordinates (default: subgraphCoordinates.csv)',
                        default='subgraphCoordinates.xyz'
    )
    
    parser.add_argument(
                        '-af', '--anisotropyFactor',
                        help='Calculate anisotropy factor (default: false)',
                        action='store_true'
    )
    
    parser.add_argument(
                        '--anisotropy-file',
                        help='Output file name for anisotropy factor (default: anisotropyFactor.csv)',
                        default='anisotropyFactor.csv'
    )
    
    parser.add_argument(
                        '-as', '--asphericityParameter',
                        help='Calculate asphericity parameter (default: false)',
                        action='store_true'
    )
    
    parser.add_argument(
                        '--asphericity-file',
                        help='Output file name for asphericity parameter (default: asphericityParameter.csv)',
                        default='asphericityParameter.csv'
    )

    
    args = vars(parser.parse_args())
    
    # Check if the required arguments are provided
    if len(sys.argv) == 1:
        parser.print_help()
        print("\nNo arguments provided. Try again ...")
        sys.exit(1)
    
    return args

def element_distance_tuple(value):
    try:
        element, distance = value.split(':')
        return (element, float(distance))
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"Invalid format: '{value}'. Expected format: ELEMENT:DISTANCE (e.g., O:3.5)"
        )
          

def wrapXYZ(xyz: pd.DataFrame, x: float, y: float, z: float) -> pd.DataFrame:
    
    """
    Apply periodic boundary conditions to the XYZ coordinates.
    
    Args:
        xyz (pd.DataFrame): DataFrame containing the coordinates of atoms.
        x (float): Box length in x-direction.
        y (float): Box length in y-direction.
        z (float): Box length in z-direction.
        
    Returns:
        pd.DataFrame: DataFrame with wrapped coordinates.
    """
    
    # Wrap coordinates
    wrapped_xyz = xyz.copy()
    wrapped_xyz['x'] = wrapped_xyz['x'] % x
    wrapped_xyz['y'] = wrapped_xyz['y'] % y
    wrapped_xyz['z'] = wrapped_xyz['z'] % z
    
    return wrapped_xyz