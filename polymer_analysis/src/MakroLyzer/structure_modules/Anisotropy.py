import numpy as np

from MakroLyzer.structure_modules.structureBase import StructureAnalyzer
from MakroLyzer.structure_modules.Helper import get_Gtensor_eigVal_eigVec

class AnisotropyAnalyzer(StructureAnalyzer):
    """
    Analyzer for the calculation of the Anisotropy of a given molecular graph.
    Inherits from StructureAnalyzer.
    
    The anisotropy factor is calculated as follows:
    (Arkin, H.; Janke, W. Ground-state properties of a polymer chain in an attractive sphere. J. Phys. Chem. B 2012, 116, 10379-10386.)
    
    κ² = 1 - 3*(λ1λ2 + λ2λ3 + λ3λ1) / (λ1 + λ2 + λ3)²
    
    where λ1 >= λ2 >= λ3 are the eigenvalues of the radius of gyration tensor.
    The κ² values range from 0 to 1. A hexagon and a dodecahedron, for instance, are highly symmetrical two and three-dimensional objects, 
    with κ² values of 0.25 and 0, respectively. Linear chains, on the other hand, exhibit a κ² = 1. 
    """
    
    def __init__(self, output_handler=None):
        """
        Initialize AnisotropyAnalyzer.

        Args:
            output_handler (OutputHandler): Handler for writing output.
        """
        super().__init__(output_handler)
        
    def initialize_output(self):
        """
        Initialize output file with header. ('streaming' mode)
        """
        if self.output_handler:
            header = "Frame, Anisotropy Factor (κ²)"
            if self.output_handler.mode == 'streaming':
                self.output_handler.initialize_file(header)
                
    def compute(self, graph):
        """
        Calculate the Anisotropy Factor for the given graph.

        Args:
            graph (GraphManager):Graph to analyze.
        """
        # get eigenvalues and sort in descending order
        eigVal, _ = get_Gtensor_eigVal_eigVec(graph)
        eigVal = np.sort(eigVal)[::-1]
        
        kappa_sq = 1 - 3 * (eigVal[0] * eigVal[1] + eigVal[1] * eigVal[2] + eigVal[2] * eigVal[0]) / (eigVal.sum() ** 2)
        
        return kappa_sq
    
    def render_output(self, data, frame_idx):
        """
        Write/Save data for this frame.
        
        Args:
            data (list): The computed kappa_sq data.
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
        header = "Frame, Anisotropy Factor (κ²)"
        super().finalize_output(header)
        
def get_anisotropy_factor(graph):
    """
    Calculate the Anisotropy Factor for the given graph.

    Args:
        graph (GraphManager): Graph to analyze.
        
    Returns:
        float: The computed anisotropy factor (κ²).
    """
    analyzer = AnisotropyAnalyzer()
    kappa_sq = analyzer.compute(graph)
    return kappa_sq
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        