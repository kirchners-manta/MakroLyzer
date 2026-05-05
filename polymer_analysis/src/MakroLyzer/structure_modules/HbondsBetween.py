from MakroLyzer.structure_modules.structureBase import StructureAnalyzer
from MakroLyzer.input_handling import subgraphSelection

class HbondsBetweenAnalyzer(StructureAnalyzer):
    """
    Analyzer for hydrogen bonds in a graph between two selections.
    Inherits from StructureAnalyzer.
    Hanldes computation and output of HBond counts per frame.
    """

    requires_full_graph = True 

    def __init__(self, cutoffs, selection1, selection2, output_handler=None, static_topology=False):
        """
        Initialize HbondsBetween analyzer.

        Args:
            cutoffs (list): List of tuples (element type, max. H-acceptor dist,
                max donor-acceptor dist, angle cutoff).
            selection1 (str): First subgraph selection token.
            selection2 (str): Second subgraph selection token.
            output_handler (OutputHandler): Handler for writing output.
            static_topology (bool): Cache selected nodes after the first frame.
        """
        super().__init__(output_handler)
        self.cutoffs = cutoffs
        self.selection1_label = selection1
        self.selection2_label = selection2
        self.selection1 = subgraphSelection.parse_selection_list([selection1])
        self.selection2 = subgraphSelection.parse_selection_list([selection2])
        self.static_topology = static_topology
        self.selection1_nodes = None
        self.selection2_nodes = None

    def initialize_output(self):
        """
        Initialize output file with header. ('streaming' mode)
        """
        if self.output_handler:
            header = (
                "Frame, Selection 1, Selection 2, Element Type, H-Acceptor dist / Å, "
                "Donor-Acceptor dist / Å, Angle cutoff / °, Number of Hydrogen Bonds"
            )
            if self.output_handler.mode == 'streaming':
                self.output_handler.initialize_file(header)

    def _get_selection_nodes(self, graph):
        if self.selection1_nodes is None or self.selection2_nodes is None or not self.static_topology:
            self.selection1_nodes = graph.select_subgraph_nodes(self.selection1)
            self.selection2_nodes = graph.select_subgraph_nodes(self.selection2)
        return self.selection1_nodes, self.selection2_nodes

    def compute(self, graph):
        """
        Calculate the number of H-bonds between the two selected species.

        Args:
            graph (GraphManager): Full graph to analyze.

        Returns:
            List of tuples (selection1, selection2, elementType, HAcceptor_dist,
            DonorAcceptor_dist, Angle_cut, count).
        """
        selection1_nodes, selection2_nodes = self._get_selection_nodes(graph)
        numberOfHbonds = []

        for elementType, HAcceptor_dist, DonorAcceptor_dist, Angle_cut in self.cutoffs:
            hbonds = graph.get_hbonds_between(
                elementType,
                HAcceptor_dist,
                DonorAcceptor_dist,
                Angle_cut,
                selection1_nodes,
                selection2_nodes,
            )
            numberOfHbonds.append(len(hbonds))

        return [
            (
                self.selection1_label,
                self.selection2_label,
                elementType,
                HAcceptor_dist,
                DonorAcceptor_dist,
                Angle_cut,
                num,
            )
            for (elementType, HAcceptor_dist, DonorAcceptor_dist, Angle_cut), num in zip(self.cutoffs, numberOfHbonds)
        ]

    def render_output(self, data, frame_idx):
        """
        Write/Save H-bond-between data for this frame.

        Args:
            data (list): The computed H-bond-between data.
            frame_idx (int): Current frame number.
        """
        for selection1, selection2, element_type, HAcceptor_dist, DonorAcceptor_dist, Angle_cut, count in data:
            row = (
                f"{frame_idx},"
                f"{selection1},"
                f"{selection2},"
                f"{element_type},"
                f"{HAcceptor_dist:.3f},"
                f"{DonorAcceptor_dist:.3f},"
                f"{Angle_cut:.3f},"
                f"{count}"
            )
            self.output_handler.append_row(row)

    def finalize_output(self):
        """
        Finalize output file (write header and rows - 'collect' mode).
        """
        header = (
            "Frame, Selection 1, Selection 2, Element Type, H-Acceptor dist / Å, "
            "Donor-Acceptor dist / Å, Angle cutoff / °, Number of Hydrogen Bonds"
        )
        super().finalize_output(header)
