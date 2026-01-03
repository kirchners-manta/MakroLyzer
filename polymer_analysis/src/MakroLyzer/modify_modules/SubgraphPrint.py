from MakroLyzer.modify_modules.structureModifierBase import StructureModifier

class SubgraphPrintModifier(StructureModifier):
    """
    Modifier to print subgraphs of molecular graphs.
    Inherits from StructureModifier.
    """

    def __init__(self, output_handler=None):
        """
        Initialize SubgraphPrintModifier.

        Args:
            output_handler (OutputHandler): Handler for writing output.
        """
        super().__init__(output_handler)

    def compute(self, graph):
        """
        Get subgraphs of the molecular graph.

        Args:
            graph : The molecular graph to modify.

        Returns:
            Subgraphs of the molecular graph.
        """
        subgraphs = list(graph.get_subgraphs())
        return subgraphs

    def render_output(self, subgraphs, frame_idx):
        """
        Write graph xyz and fragment IDs to file.

        Args:
            subgraphs (list): List of subgraphs.
            frame_idx (int): Current frame number.
        """
        for i, subgraph in enumerate(subgraphs):
            i += 1  # Start fragment IDs at 1
            NoNodes = subgraph.number_of_nodes()
            header = f"{NoNodes}\nelement  x         y        z"
            self.output_handler.write_output(header, subgraph, ID=i)
