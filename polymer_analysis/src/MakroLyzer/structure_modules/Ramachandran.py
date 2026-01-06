from MakroLyzer.structure_modules.structureBase import StructureAnalyzer

class RamachandranAnalyzer(StructureAnalyzer):
    """
    Analyzer for the calculation of Ramachandran matricies. 
    Inherits from StructureAnalyzer.
    
    A Ramachandran matrix contains the phi/psi counts of the backbone of an amino acid.
    
    phi : C(i)-N(i+1)-Ca(i+1)-C(i+1)
    psi : N(i)-Ca(i)-C(i+1)-N(i+1)
    
    The backbone of an amino acid is build like this:
    ...-N-Calpha-CarbonylC-N-...
    """
    
    def __init__(self, output_handler=None):
        """
        Initialize the RamachandranAnalyzer.

        Args:
            output_handler (OutputHandler): Handler for writing output.
        """
        super().__init__(output_handler)
        
    def initialize_output(self):
        pass
    
    def compute(self, graph):
        """
        Get the Ramachandran matrix of the graph.
        The Polymer needs to be an Amino Acid or an alternative Amino Acid. 

        Args:
            graph (GraphManager): The graph to calculate the Ramachandran plot data for.

        Returns:
            Matrix (list of lists): A Ramachandran matrix containing the phi and psi angle counts.
        """  

        # Prepare a ramachandran matrix of ints and initialize to 0
        ramachandran = [[0 for _ in range(360)] for _ in range(360)]

        # Get subgraphs
        subgraphs = graph.get_subgraphs()

        for subgraph in subgraphs:

            # Get the backbone of the AA which is always like this: 
            # N-C-C-N-...
            # N-Calpha-CarbonylC-N-...
            backbone = subgraph.AminoAcidBackbone()

            # Slide over the backbone and calculate phi/psi angles
            # Phi is defined as the dihedral angle C(i)-N(i+1)-Ca(i+1)-C(i+1)
            # Psi is defined as the dihedral angle N(i)-Ca(i)-C(i+1)-N(i+1)

            # We start with the first dihedral at the third atom in the backbone which is the
            # first CarbonylC = Cprime
            i = 2
            while i < len(backbone) - 4:
                Cprime1 = backbone[i]
                N1 = backbone[i+1]
                Calpha = backbone[i+2]
                Cprime2 = backbone[i+3]
                N2 = backbone[i+4]

                # phi C(i)-N(i+1)-Ca(i+1)-C(i+1)
                phi = graph.dihedral(Cprime1, N1, Calpha, Cprime2, sign=True)
                phi = round(phi) + 180
                # psi N(i)-Ca(i)-C(i+1)-N(i+1)
                psi = graph.dihedral(N1, Calpha, Cprime2, N2, sign=True)
                psi = round(psi) + 180

                # Increment the ramachandran matrix
                ramachandran[phi][psi] += 1

                # move to the next Calpha
                i += 3 
                
        return ramachandran
       
    def render_output(self, data, frame_idx):
        """
        Write/Save data for this frame.
        
        Args:
            data (list): The computed b data.
            frame_idx (int): Current frame number.
            """
        self.output_handler.append_matrix(data, frame_idx)
        
    def finalize_output(self):
        """
        Finalize output after all frames are processed ('collect' mode)
        """
        if self.output_handler and self.output_handler.mode == 'collect':
            self.output_handler.finalize_matrices()