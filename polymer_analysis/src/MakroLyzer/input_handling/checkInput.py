import os

class FileNotFoundError(Exception):
    """Exception raised when a file is not found."""
    pass

class InvalidFileFormatError(Exception):
    """Exception raised when a file format is invalid."""
    pass

class EmptyFileError(Exception):
    """Exception raised when a file is empty."""
    pass

def checkInput(args):
    """
    Check the input arguments and validate the file paths.
    
    Args:
        args (dict): Dictionary containing the input arguments.
        
    Raises:
        FileNotFoundError: If the specified file is not found.
        InvalidFileFormatError: If the file format is invalid.
    """
    
    #---XYZ---#
    # Check if the XYZ file exists
    if args['xyzFile'] is not None:
        xyzFilePath = args['xyzFile']
        if not os.path.isfile(xyzFilePath):
            raise FileNotFoundError(f"XYZ file '{xyzFilePath}' not found.")

        # Check if the file format is valid
        if not xyzFilePath.endswith('.xyz'):
            raise InvalidFileFormatError(f"Invalid file format for '{xyzFilePath}'. Expected .xyz file.")

        # Check if the file is empty
        if os.path.getsize(xyzFilePath) == 0:
            raise EmptyFileError(f"File '{xyzFilePath}' is empty.")
        
    #---LAMMPS---#
    if args['lmpFile'] is not None:
        lmpFilePath = args['lmpFile']
        if not os.path.isfile(lmpFilePath):
            raise FileNotFoundError(f"LAMMPS file '{lmpFilePath}' not found.")

        # Check if the file format is valid
        if not lmpFilePath.endswith('.lmp') and not lmpFilePath.endswith('.lammpstrj'):
            raise InvalidFileFormatError(f"Invalid file format for '{lmpFilePath}'. Expected .lmp or .lammpstrj file.")

        # Check if the file is empty
        if os.path.getsize(lmpFilePath) == 0:
            raise EmptyFileError(f"File '{lmpFilePath}' is empty.")
        
    
    # GROMACS trajectories need atom identities from a separate topology.
    binary_trajectory = args.get('trrFile') or args.get('xtcFile')
    if binary_trajectory and not args.get('topologyFile'):
        raise ValueError('TRR/XTC input requires -tpr/--topology (TPR, GRO, or PDB).')
    if args.get('topologyFile') and not binary_trajectory:
        raise ValueError('-tpr/--topology requires a TRR or XTC trajectory.')
    for key, extensions in (
        ('groFile', ('.gro',)), ('pdbFile', ('.pdb',)),
        ('trrFile', ('.trr',)), ('xtcFile', ('.xtc',)),
        ('topologyFile', ('.tpr', '.gro', '.pdb')),
    ):
        path = args.get(key)
        if path is None:
            continue
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Input file '{path}' not found.")
        if not str(path).lower().endswith(extensions):
            raise InvalidFileFormatError(
                f"Invalid file format for '{path}'. Expected {' or '.join(extensions)} file."
            )
        if os.path.getsize(path) == 0:
            raise EmptyFileError(f"File '{path}' is empty.")

    #---Repeating Units---#
    if args['patternFile'] is not None:
        # Check if the pattern file exists
        patternPath = args['patternFile']
        if not os.path.isfile(patternPath):
            raise FileNotFoundError(f"Pattern file '{patternPath}' not found.")

        # Check if the file format is valid
        if not patternPath.endswith('.txt'):
            raise InvalidFileFormatError(f"Invalid file format for '{patternPath}'. Expected .txt file.")

        # Check if the file is empty
        if os.path.getsize(patternPath) == 0:
            raise EmptyFileError(f"File '{patternPath}' is empty.")

    #---Wrap Graph Print---#
    if args.get('wrapGraphPrint') is not None and args.get('BoxSize') is None:
        raise ValueError("'-wrap/--wrapGraphPrint' requires '-bs/--BoxSize'.")
        
    
