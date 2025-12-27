import numpy as np

from MakroLyzer.structure_modules.structureBase import StructureAnalyzer
from MakroLyzer.structure_modules.Helper import get_Gtensor_eigVal_eigVec

class AsphericityAnalyzer(StructureAnalyzer):
    """
    Analyzer for the calculation of the Asphericity of a given molecular graph.
    Inherits from StructureAnalyzer. 
    
    The asphericity parameter is calculated as follows:
    (Arkin, H.; Janke, W. Ground-state properties of a polymer chain in an attractive sphere. J. Phys. Chem. B 2012, 116, 10379-10386.)
     
        b = λ1 - 1/2 * (λ2 + λ3)
        
    where λ1 >= λ2 >= λ3 are the eigenvalues of the radius of gyration tensor.
    Values close to 0 indicate a high sphericity, while large values indicate a low sphericity.
    """
    
    def __init__(self, output_handler=None):
        """
        Initialize the AsphericityAnalyzer.

        Args:
            output_handler (OutputHandler): Handler for writing output.
        """
        super().__init__(output_handler)
        
    def initialize_output(self):
        """
        Initialize output file with header. ('streaming' mode)
        """
        if self.output_handler:
            header = "Frame, Asphericity Parameter (b)"
            if self.output_handler.mode == 'streaming':
                self.output_handler.initialize_file(header)
                
    def compute(self, graph):
        """
        Calculate the Asphericity Parameter for the given graph.

        Args:
            graph (GraphManager): Graph to analyze.
        """
        
        # get eigenvalues and sort in descending order
        eigVal, _ = get_Gtensor_eigVal_eigVec(graph)
        eigVal = np.sort(eigVal)[::-1]
        
        b = eigVal[0] - 0.5 * (eigVal[1] + eigVal[2])
        
        return b
    
    def render_output(self, data, frame_idx):
        """
        Write/Save data for this frame.
        
        Args:
            data (list): The computed b data.
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
        header = "Frame, Asphericity Parameter (b)"
        super().finalize_output(header)
        
def get_anisotropy_factor(graph):
    """
    Calculate the Anisotropy Factor for the given graph.

    Args:
        graph (GraphManager): Graph to analyze.
        
    Returns:
        float: The computed anisotropy factor (κ²).
    """
    analyzer = AsphericityAnalyzer()
    kappa_sq = analyzer.compute(graph)
    return b
        