from collections import deque

import numpy as np

from MakroLyzer.dynamic_modules.dynamicBase import DynamicAnalyzer


class RMSDAnalyzer(DynamicAnalyzer):
    """
    Analyzer for the calculation of atomic root mean square deviation averaged 
    over atoms and time origins.
    Inherits from DynamicAnalyzer.

    RMSD= sqrt(1/n sum |x_i-x_i^{ref}|^2)
    ref is the first structure
    """

    header = "Correlation Time, RMSD / Å"

    def __init__(self, output_handler=None):
        """
        Initialize the RMSDAnalyzer.
        
        Args:
            output_handler (OutputHandler): Handler for writing output.
        
        Set maximum lag time and interval between analyzed frames.

        Both times use the same units. Only positive sampled lags no greater
        than correlation_depth are included; unsampled lags are omitted.
        """
        super().__init__(output_handler)

        self.refstruc=None
        self._nodes = None

    def set_first_structure(self, coords):
        """
        Saves the first structure as reference points
        """
        self.refstruc=coords.copy()

    def initialize_output(self):
        """
        Initialize output file with header. ('streaming' mode)
        """
        if self.output_handler is not None and self.output_handler.mode == 'streaming':
            self.output_handler.initialize_file(self.header)

    def compute(self, graph):
        """
        Accumulate one atom-averaged sample per available lag.
        """
        nodes, coords = graph.get_all_coordinates()
        nodes = tuple(nodes)
        coords = np.asarray(coords, dtype=float)
        if not nodes or coords.shape != (len(nodes), 3) or not np.all(np.isfinite(coords)):
            raise ValueError("RMSD requires a nonempty set of finite 3D coordinates.")
        if self._nodes is not None and nodes != self._nodes:
            raise ValueError("RMSD atom identities and ordering must remain unchanged.")
        if self.refstruc is None:
            self.set_first_structure(coords)
        self._nodes = nodes
        self.results = None
        rmsd = np.sqrt(np.mean(np.sum((coords - self.refstruc) ** 2, axis=1)))

        return rmsd

    def render_output(self, data, frame_idx):
        """
        Write/Save data for this frame.
        
        Args:
            data (list): The computed rmsd data.
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
        header = "Frame, Root mean square deviation (RMSD)"
        super().finalize_output(header)
