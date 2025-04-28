import pandas as pd
import ast

def readPattern(patternFilePath: str) -> pd.DataFrame:
    """
    Read pattern file and extract repeating units (monomers).
    
    Args:
        patternFilePath (str): Path to the pattern file.
        
    Returns:
        pd.DataFrame: DataFrame containing the repeating units (monomers)
        and optionally the element from which side to start.
    """    
    
    data = []
    with open(patternFilePath) as file:
        lines = [line.strip() for line in file if line.strip()]  # Ignore empty lines
        
        if len(lines) == 0:
            raise ValueError("The pattern file is empty.")
        if len(lines) > 2:
            raise ValueError("The pattern file contains more than two lines.")
        
        # Safely evaluate the first line to a Python list
        try:
            pattern = ast.literal_eval(lines[0])
            if not (isinstance(pattern, list) and all(isinstance(p, list) for p in pattern)):
                raise ValueError("Pattern must be a list of lists.")
        except (SyntaxError, ValueError) as e:
            raise ValueError(f"Could not parse pattern line: {e}")
        
        # Second line: element (optional)
        element = lines[1] if len(lines) == 2 else None
    
    data.append({'pattern': pattern, 'element': element})
    df = pd.DataFrame(data)
    return df
