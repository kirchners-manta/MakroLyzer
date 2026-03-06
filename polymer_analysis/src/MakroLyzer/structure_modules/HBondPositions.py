from pathlib import Path

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
        # Per-frame files are written in render_output.
        return None

    @staticmethod
    def _midpoint(h_pos, acceptor_pos):
        h_pos = np.asarray(h_pos, dtype=float)
        acceptor_pos = np.asarray(acceptor_pos, dtype=float)
        return 0.5 * (h_pos + acceptor_pos)

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
        Write one XYZ file for this frame, containing H-bond midpoint coordinates.
        """
        if not self.output_handler:
            return

        out_path = Path(self.output_handler.file_path)
        frame_path = out_path.with_name(f"{out_path.stem}_frame{frame_idx}{out_path.suffix or '.xyz'}")

        with open(frame_path, "w") as fh:
            fh.write(f"{len(data)}\n")
            fh.write(f"HBond midpoints frame {frame_idx}\n")
            for point in data:
                x, y, z = np.asarray(point, dtype=float)
                fh.write(f"Hb {x:.6f} {y:.6f} {z:.6f}\n")

    def finalize_output(self):
        # Nothing to finalize: files are written per frame.
        return None


# Backward-compatible alias (legacy typo / naming overlap).
HBondsAnalyzer = HBondPositionsAnalyzer
