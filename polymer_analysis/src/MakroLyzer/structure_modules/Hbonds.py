from MakroLyzer.structure_modules.structureBase import StructureAnalyzer

class HBondsAnalyzer(StructureAnalyzer):
    """
    Analyzer for hydrogen bonds in a polymer graph.
    Inherits from StructureAnalyzer.
    Hanldes computation and output of HBond counts per frame.
    """
    
    def __init__(self, cutoffs, output_handler=None):
        """
        Initialize Hbonds analyzer.

        Args:
            cutoffs (list): List of tuples (element type, max. H-acceptor dist, max donor-acceptor dist, angle cutoff)
            output_handler (OutputHandler): Handler for writing output.
        """
        
        # Initialize OutputHandler
        super().__init__(output_handler)
        self.cutoffs = cutoffs
        
    def initialize_output(self):
        """
        Initialize output file with header. ('streaming' mode)
        """
        if self.output_handler:
            header = "Frame, Element Type, H-Acceptor dist / Å, Donor-Acceptor dist / Å, Angle cutoff / °, Number of Hydrogen Bonds"
            if self.output_handler.mode == 'streaming':
                self.output_handler.initialize_file(header)
                
    def compute(self, graph):
        """
        Calculate the number of Hydrogen Bonds in the given molecular graph. 

        Args:
            graph (GraphManager): Graph to analyze.
            
        Returns: 
            List of tuples (elementType, HAcceptor_dist, DonorAcceptor_dist, Angle_cut, count)
        """
        
        numberOfHbonds = []

        for elementType, HAcceptor_dist, DonorAcceptor_dist, Angle_cut in self.cutoffs:        
            # Get the hydrogen bonds for the element type
            hbonds = graph.get_hbonds(elementType, HAcceptor_dist, DonorAcceptor_dist, Angle_cut)
            numberOfHbonds.append(len(hbonds))

        hbonds = [
            (elementType, HAcceptor_dist, DonorAcceptor_dist, Angle_cut, num) 
                  for (elementType, HAcceptor_dist, DonorAcceptor_dist, Angle_cut), num in zip(self.cutoffs, numberOfHbonds)
        ]
        return hbonds
    
    def render_output(self, data, frame_idx):
        """
        Write/Save Hbond data for this frame.
        
        Args:
            data (list): The computed Hbond data.
            frame_idx (int): Current frame number.
        """
        for element_type, HAcceptor_dist, DonorAcceptor_dist, Angle_cut, count in data:
            row = (
                f"{frame_idx},"
                f"{element_type},"
                f"{HAcceptor_dist:.3f},"
                f"{DonorAcceptor_dist:.3f},"
                f"{Angle_cut:.3f},"
                f"{count}"
            )
            self.output_handler.append_row(row)
            
    def finalize_output(self):
        """
        Finalize output file (write header and rows - 'collect' mode)
        """
        # Run the parent class
        header = "Frame, Element Type, H-Acceptor dist / Å, Donor-Acceptor dist / Å, Angle cutoff / °, Number of Hydrogen Bonds"
        super().finalize_output(header)
        
def get_Hbonds(graph, cutoffs):
    """
    Calculate Hbonds for the given graph.

    Args:
        graph (GraphManager): The graph to calculate the Hbond count for.
        cutoffs (list): List of tuples (element type, max. H-acceptor dist, max donor-acceptor dist, angle cutoff)
        
    Returns: 
        list: A list of tuples containing numbers of Hbonds and their atom pairs.
    """
    analyzer = HBondsAnalyzer(cutoffs)
    return analyzer.compute(graph)