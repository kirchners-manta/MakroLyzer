import numpy as np

from MakroLyzer.structure_modules.structureBase import StructureAnalyzer


class HBondPositionsAnalyzer(StructureAnalyzer):
    """
    Analyzer that writes one XYZ file per frame with H-bond midpoint positions.
    """

    def __init__(self, cutoffs, output_handler=None):
        super().__init__(output_handler)
        self.cutoffs = cutoffs

    def initialize_output(self):
        """
        No global output file to initialize.
        One XYZ file is generated per frame.
        """
        pass

    @staticmethod
    def _midpoint(h_pos, acceptor_pos):
        h_pos = np.asarray(h_pos, dtype=float)
        acceptor_pos = np.asarray(acceptor_pos, dtype=float)
        return 0.5 * (h_pos + acceptor_pos)

    @staticmethod
    def _to_xyz_matrix(points, frame_idx):
        """
        Convert midpoint coordinates to XYZ file rows for one frame.
        """
        matrix = [[str(len(points))], [f"HBond midpoints frame {frame_idx}"]]
        for point in points:
            x, y, z = np.asarray(point, dtype=float)
            matrix.append([f"Hb {x:.6f} {y:.6f} {z:.6f}"])
        return matrix

    def compute(self, graph):
        """
        Returns:
            list[np.ndarray]: Midpoint coordinates for all H-bonds in this frame.
        """
        points = []
        for element_type, h_acceptor_dist, donor_acceptor_dist, angle_cut in self.cutoffs:
            hbonds = graph.get_hbonds(
                element_type,
                h_acceptor_dist,
                donor_acceptor_dist,
                angle_cut,
            )
            for node_h, node_acceptor in hbonds:
                h_pos = graph.get_coordinates(node_h)
                acceptor_pos = graph.get_coordinates(node_acceptor)
                points.append(self._midpoint(h_pos, acceptor_pos))
        return points

    def render_output(self, data, frame_idx):
        """
        Write/save one XYZ file for this frame.
        """
        matrix = self._to_xyz_matrix(data, frame_idx)
        if not self.output_handler:
            return

        if self.output_handler.mode == "streaming":
            self.output_handler.append_matrix(matrix, frame_idx)
            return

        # Keep per-frame file behavior consistent with historical implementation:
        # write each frame file immediately even when handler mode is "collect".
        original_mode = self.output_handler.mode
        self.output_handler.mode = "streaming"
        try:
            self.output_handler.append_matrix(matrix, frame_idx)
        finally:
            self.output_handler.mode = original_mode

    def finalize_output(self):
        """
        Finalize output after all frames are processed.
        """
        return None
