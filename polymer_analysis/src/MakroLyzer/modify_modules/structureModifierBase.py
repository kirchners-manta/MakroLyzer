"""
This is the base class for all structure modifier modules.
Each structure modifier handles computation and its own output formatting.

-> All structure modules inherit from this class. 
   Thus, they implement the methods defined here.
"""

# abc - Abstract Base Classes
# Infrastructure for defining abstract base classes in Python
from abc import ABC, abstractmethod
# os - Miscellaneous operating system interfaces
# Provides a way of using operating system dependent functionality
import os


class ModifyOutputHandler:
    """
    Base class for handling the output of a structure modifier.
    """
    
    def __init__(self, file_path):
        """
        Initialize Output Handler.

        Args:
            file_path (str): Path where to save the output file.
        """
        self.file_path = file_path
        
    def write_output(self, header, graph, attributes=None):
        """
        Write fragment data for `graph` to CSV.

        Args:
            header (str): Header row for the CSV file.
            graph: GraphManager instance whose node data will be written. (xyz coordinates)
            attributes (iterable, optional): Additional node attributes to include as extra columns.
        """
        # If the file exists, append a number to the filename to avoid overwriting
        base, ext = os.path.splitext(self.file_path)
        counter = 1
        while os.path.exists(self.file_path):
            self.file_path = f"{base}_{counter}{ext}"
            counter += 1
                   
        # Write data           
        with open(self.file_path, 'w') as f:
            f.write(header + '\n')
            for node in graph.nodes():
                d = graph.nodes[node]
                elem = d.get('element', 'X')
                x = d.get('x', 0.0)
                y = d.get('y', 0.0)
                z = d.get('z', 0.0)
                row = f"{elem:6}{x:10.4f}{y:10.4f}{z:10.4f}"
                if attributes:
                    for attr in attributes:
                        row += f"    {d.get(attr, '')}"
                f.write(row + '\n')                                
                    
class StructureModifier(ABC):
    """
    Base class for all structure modifiers.
    Inherits from ABC to enforce implementation of abstract methods.
    Subclasses must implement compute() and render_output() methods.
    """
    
    def __init__(self, output_handler = None):
        """
        Initialize structure modifier.

        Args:
            output_handler (OutputHandler): Handler for writing output. 
        """
        self.output_handler = output_handler
        self.frame_number = 0
        
    @abstractmethod
    def compute(self, graph):
        """
        Perform some creazy calculations to obtain the results. 

        Args:
            graph : The molecular graph to modify.
            
        Returns; 
            Results (format depends on the modifier.)
        """
        pass
    
    @abstractmethod
    def render_output(self, data, frame_idx):
        """
        Format and write the output for a single frame.

        Args:
            data : Results from the compute function.
            frame_idx (int): Current frame number.
            
        Returns:
            Formatted output.
        """
        pass
            
    def run(self, graph, frame_idx):
        """
        Execute modification for a single frame: compute + output

        Args:
            graph : The molecular graph to modify
            frame_idx (int): Current frame number
            
        Returns: 
           Result.
        """
        result = self.compute(graph)
        if self.output_handler:
            self.render_output(result, frame_idx)
        return result