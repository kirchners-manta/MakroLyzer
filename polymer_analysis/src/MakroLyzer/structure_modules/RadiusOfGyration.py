import numpy as np

from MakroLyzer.structure_modules.structureBase import StructureAnalyzer

class RadiusOfGyrationAnalyzer(StructureAnalyzer):
    """
    Analyzer for the Radius of Gyration (Rg) of a given molecular graph.
    Inherits from StructureAnalyzer.
    
    The Radius of Gyration is calculated as the square root of the average of the 
    squared distances of each atom from the center of mass.
    
    Rg² = 1/N * Σ_(i=1)^(N) (r_i - r_com)²
    
    where N is the total number of atoms, r_i is the position vector of atom i and 
    r_com is the position vector of the center of mass.
    """
    
    def __init__(self, output_handler=None):
        """
        Initialize RadiusOfGyrationAnalyzer.

        Args:
            output_handler (OutputHandler): Handler for writing output.
        """
        super().__init__(output_handler)
        
    def initialize_output(self):
        """
        Initialize output file with header. ('streaming' mode)
        """
        if self.output_handler:
            header = "Frame, Rg / Å"
            if self.output_handler.mode == 'streaming':
                self.output_handler.initialize_file(header)
                
    def compute(self, graph):
        """
        Calculate the Rg for the given graph.

        Args:
            graph (GraphManager):Graph to analyze.
        """
        
        _, coords = graph.get_all_coordinates()
        com = graph.get_com()
        squared_distances = np.sum((coords-com)**2, axis=1) # axis=1 -> sum over rows
        Rg = np.sqrt(np.mean(squared_distances))
        
        return Rg
    
    def render_output(self, data, frame_idx):
        """
        Write/Save data for this frame.
        
        Args:
            data (list): The computed Rg data.
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
        header = "Frame, Rg / Å"
        super().finalize_output(header)
        
def get_radius_of_gyration(graph):
    """
    Calculate the Radius of Gyration for the given graph.

    Args:
        graph (GraphManager): Graph to analyze.
        
    Returns:
        float: The computed Radius of Gyration.
    """
    analyzer = RadiusOfGyrationAnalyzer()
    Rg = analyzer.compute(graph)
    return Rg