import numpy as np

from MakroLyzer.structure_modules.structureBase import StructureAnalyzer

class SpreadRadiusAnalyzer(StructureAnalyzer):
    """
    Analyzer for the Spread Radius (Rs) of a given molecular graph.
    Inherits from StructureAnalyzer.
    
    The Spread Radius is calculated as the distance of the most distant atom
    from the center of mass.
    """
    
    def __init__(self, output_handler=None):
        """
        Initialize SpreadRadiusAnalyzer.

        Args:
            output_handler (OutputHandler): Handler for writing output.
        """
        super().__init__(output_handler)
        
    def initialize_output(self):
        """
        Initialize output file with header. ('streaming' mode)
        """
        if self.output_handler:
            header = "Frame, Rs / Å"
            if self.output_handler.mode == 'streaming':
                self.output_handler.initialize_file(header)
                
    def compute(self, graph):
        """
        Calculate the Rs for the given graph.

        Args:
            graph (GraphManager):Graph to analyze.
        """
        
        _, coords = graph.get_all_coordinates()
        com = graph.get_com()
        Rs = 0
        for position in coords:
            dist = np.linalg.norm(position - com)
            if dist > Rs:
                Rs = dist
            
        return Rs
    
    def render_output(self, data, frame_idx):
        """
        Write/Save data for this frame.
        
        Args:
            data (list): The computed Rs data.
            frame_idx (int): Current frame number.
        """
        row = (
            f"{frame_idx},"
            f"{data:.3f}"
        )
        self.output_handler.append_row(row)
        
    def finalize_output(self):
        """
        Finalize output file (write header and rows - 'collect' mode)
        """
        header = "Frame, Rs / Å"
        super().finalize_output(header)
        