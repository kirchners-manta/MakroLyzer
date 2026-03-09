from scipy.spatial import ConvexHull, QhullError

from MakroLyzer.structure_modules.structureBase import StructureAnalyzer

class ConvexHullVolumeAnalyzer(StructureAnalyzer):
    """
    Analyzer for the calculation of the volume of the Convex Hull spanned by a 
    particle. 
    Inherits from StructureAnalyzer.
    """
    
    def __init__(self, output_handler=None):
        """
        Initialize ConvexHullVolumeAnalyzer.

        Args:
             output_handler (OutputHandler): Handler for writing output.
        """
        super().__init__(output_handler)
        
    def initialize_output(self):
        """
        Initialize output file with header. ('streaming' mode)
        """
        if self.output_handler:
            header = "Frame, Convex Hull - Volume / Å³, Mass / g/mol, Density / g/cm³"
            if self.output_handler.mode == 'streaming':
                self.output_handler.initialize_file(header)
            
    def compute(self, graph):
        """
        Calculate the convex hull and obtain the volume within it.
        
        Args:
            graph (GraphManager):Graph to analyze.
        """        
        _, coords = graph.get_all_coordinates()

        # Convex hull volume in 3D needs at least 4 non-coplanar points. Angstrom³ 
        if len(coords) < 4:
            volume = 0.0
        else:
            try:
                hull = ConvexHull(coords)
                volume = float(hull.volume)
            except QhullError:
                volume = 0.0
                
        # Mass of the polymer. g/mol
        mass = graph.get_mass()
        
        # Density in g/cm³ (convert from g/mol and Å³)
        if volume > 0:
            density = (mass / volume) * 1e24 / 6.022e23  # g/cm³
        else:            
            density = 0.0

        return volume, mass, density

    def render_output(self, data, frame_idx):
        """
        Write/Save data for this frame.
        
        Args:
            data (list): The computed volume.
            frame_idx (int): Current frame number.
        """
        row = (
            f"{frame_idx},"
            f"{data[0]:.3f},"
            f"{data[1]:.3f},"
            f"{data[2]:.3f}"
        )
        self.output_handler.append_row(row)
        
    def finalize_output(self):
        """
        Finalize output file (write header and rows - 'collect' mode)
        """
        header = "Frame, Convex Hull - Volume / Å³, Mass / g/mol, Density / g/cm³"
        super().finalize_output(header)
        