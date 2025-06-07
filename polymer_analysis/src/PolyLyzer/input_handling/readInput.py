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
                        help='Path to the XYZ file', 
                        required=True) 
    
    parser.add_argument('-p', '--patternFile', 
                        help='Path to the TXT file -> Finds repeating units (default: false)')
    
    parser.add_argument('-s', '--saturation', 
                        help='Saturate the ends polymers (default: false)', 
                        action='store_true')
    
    parser.add_argument('-pbc', '--PBC_xyz', 
                        nargs=3, 
                        type=float, 
                        help='Apply periodic boundary conditions in x, y, z (e.g., 150 150 150)')
    
    parser.add_argument('-e2e', '--endToEndDistance', 
                        help='Calculate end-to-end distance (default: false)', 
                        action='store_true')
    
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
                        

    
    args = vars(parser.parse_args())
    
    # Check if the required arguments are provided
    if len(sys.argv) == 1:
        parser.print_help()
        print("\nNo arguments provided. Try again ...")
        sys.exit(1)
    
    return args


def readXYZ(xyzFilePath: str) -> pd.DataFrame:
    
    """
    Read XYZ file and extract coordinates.
    
    Args:
        xyzFilePath (str): Path to the XYZ file.
        
    Returns:
        pd.DataFrame: DataFrame containing the coordinates of atoms.
    """    
    
    data = []
    with open(xyzFilePath, 'r') as file:
        # Skip header lines
        lines = file.readlines()[2:]  
        for line in lines:
            parts = line.split()
            if len(parts) == 4: 
                atom, x, y, z = parts
                data.append([atom, float(x), float(y), float(z)])
            else:
                print(f"Invalid line format: {line.strip()}")
                continue

    # Create DataFrame from the list of data
    df = pd.DataFrame(data, columns=['atom', 'x', 'y', 'z'])
    # add a column with the index
    df['index'] = df.index
    return df            

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