
from MakroLyzer.errorOutputs.ErrorOutputs import ErrorOutputs

from MakroLyzer.structure_modules.structureBase import StructureAnalyzer

class MoleculeCountAnalyzer(StructureAnalyzer):
    """
    Analyzer for counting the Molecules in a given graph.
    Inherits from StructureAnalyzer.
    The molecules are devided in chains and rings.
    """
    
    def __init__(self, output_handler=None):
        """
        Initialize MoleculeCountAnalyzer.

        Args:
            output_handler (OutputHandler): Handler for writing the output.
        """
        super().__init__(output_handler)
        
    def initialize_output(self):
        """
        Initialize output file with header. ('streaming' mode)
        """
        if self.output_handler:
            header = "Frame, Molecule count, Chain count, Ring count"
            if self.output_handler.mode == 'streaming':
                self.output_handler.initialize_file(header)
                
    def compute(self, graph):
        """
        Calculate the number of molecules in the given molecular graph,
        as well as the chain and ring counts.

        Args:
            graph (GraphManager): Graph to analyze.
            
        Returns: 
            Number of Molecules.
        """
        NoMolecules = len(graph.get_subgraphs())
        
        NoRings = 0
        NoChains = 0
        
        # If the molecule has free ends, it is a chain, else a ring
        for molecule in graph.get_subgraphs():
            backbone = molecule.find_longest_path()
            
            # Fileter graph without 1-order nodes
            G_wo1 = molecule.remove_1order()
            backbone = G_wo1.subgraph(backbone).copy()
            
            if any(backbone.degree(n) == 1 for n in backbone.nodes()):
                NoChains += 1
            else:
                NoRings += 1
                
        total_molecules = NoRings + NoChains
        
        if total_molecules != NoMolecules:
            raise ValueError(ErrorOutputs.MOLECULE_COUNT_MISMATCH_ERROR)
        
        MoleculeCounts = [(NoMolecules, NoChains, NoRings)]
        
        return MoleculeCounts
            
    def render_output(self, data, frame_idx):
        """
        Write/Save data for this frame.
        
        Args:
            data (list): The computed molecule count.
            frame_idx (int): Current frame number.
        """
        for NoMolecules, NoChains, NoRings in data:
            row = (
                f"{frame_idx},"
                f"{NoMolecules},"
                f"{NoChains},"
                f"{NoRings}"
            )
            self.output_handler.append_row(row)
        
    def finalize_output(self):
        """
        Finalize output file (write header and rows - 'collect' mode)
        """
        header = "Frame, Molecule count, Chain count, Ring count"
        super().finalize_output(header)
        
def get_anisotropy_factor(graph):
    """
    Calculate the Anisotropy Factor for the given graph.

    Args:
        graph (GraphManager): Graph to analyze.
        
    Returns:
        float: The computed molecule count, chain and ring count.
    """
    analyzer = MoleculeCountAnalyzer()
    counts = analyzer.compute(graph)
    return counts
        
    