from MakroLyzer.modify_modules.structureModifierBase import StructureModifier

class AminoSaturationModifier(StructureModifier):
    """
    Modifier to saturate amino particles in molecular graphs.
    Inherits from StructureModifier.
    """

    def __init__(self, output_handler=None):
        """
        Initialize NylonSaturationModifier.

        Args:
            output_handler (OutputHandler): Handler for writing output.
        """
        super().__init__(output_handler)

    def compute(self, graph):
        """
        Saturate the nylon chains in the molecular graph.

        Args:
            graph : The molecular graph to modify.

        Returns:
            Results of the analysis (format depends on the modifier).
        """
        graph.saturate()
        return graph

    def render_output(self, graph, frame_idx):
        """
        Write graph xyz to file.

        Args:
            graph (GraphManager): The modified graph.
            frame_idx (int): Current frame number.
        """
        NoNodes = graph.number_of_nodes()
        header = f"{NoNodes}\nelement  x         y        z"
        self.output_handler.write_output(header, graph)        