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
    xyzFilePath = args['xyzFile']
    if not os.path.isfile(xyzFilePath):
        raise FileNotFoundError(f"XYZ file '{xyzFilePath}' not found.")
    
    # Check if the file format is valid
    if not xyzFilePath.endswith('.xyz'):
        raise InvalidFileFormatError(f"Invalid file format for '{xyzFilePath}'. Expected .xyz file.")
    
    # Check if the file is empty
    if os.path.getsize(xyzFilePath) == 0:
        raise EmptyFileError(f"File '{xyzFilePath}' is empty.")
        
    
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
        
    