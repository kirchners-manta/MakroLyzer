import os
import ast
import pandas as pd

from MakroLyzer.modify_modules.structureModifierBase import StructureModifier

class PatternModifier(StructureModifier):
    """
    Modifier assigns specific pattern IDs to subgraphs based on predefined patterns.
    Inherits from StructureModifier.
    """

    def __init__(self, pattern_file, output_handler=None):
        """
        Initialize PatternModifier.
        
        Args:
            pattern_file (str): Path to the pattern definition file.
            output_handler (OutputHandler): Handler for writing output.
        """
        super().__init__(output_handler)
        self.patterns = self._read_patterns(pattern_file)
                
    def _read_patterns(self, pattern_file):
        """
        Read predefined patterns from a file.

        Args:
            pattern_file (str): Path to the pattern definition file.
        """
        # check if file exists
        if not os.path.exists(pattern_file):
            raise FileNotFoundError(f"Pattern file {pattern_file} not found.")
        
        # Read patterns into a dictionary
        data = []
        with open(pattern_file) as file:
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
    
    def compute(self, graph):
        """ 
        Obtain IDs for each pattern in the graph.
        """
        graph.find_and_tag_patterns(
            self.patterns['pattern'].values[0],
            self.patterns['element'].values[0]
        )
        return graph
        
    def render_output(self, graph, frame_idx):
        """
        Write graph xyz and fragment IDs to file.

        Args:
            graph (GraphManager): The modified graph.
            frame_idx (int): Current frame number.
        """
        NoNodes = graph.number_of_nodes()
        header = f"{NoNodes}\nelement x          y         z          FragmentID"
        self.output_handler.write_output(header, graph, attributes=['fragmentID'])        