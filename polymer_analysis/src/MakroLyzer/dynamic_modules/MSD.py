from collections import deque

import numpy as np

from MakroLyzer.dynamic_modules.dynamicBase import DynamicAnalyzer


class MSDAnalyzer(DynamicAnalyzer):
    """
    Analyzer for the calculation of atomic mean squared displacement averaged 
    over atoms and time origins.
    Inherits from DynamicAnalyzer.

    For lag k, average mean_i(|r_i(t + k) - r_i(t)|**2) over all
    available origins t. 
    """

    header = "Correlation Time, MSD / Å²"

    def __init__(self, correlation_depth, timestep, output_handler=None):
        """
        Initialize the MSDAnalyzer.
        
        Args:
            correlation_depth (float): Maximum lag time to consider for MSD calculation.
            timestep (float): Time interval between analyzed frames.
            output_handler (OutputHandler): Handler for writing output.
        
        Set maximum lag time and interval between analyzed frames.

        Both times use the same units. Only positive sampled lags no greater
        than correlation_depth are included; unsampled lags are omitted.
        """
        super().__init__(output_handler)
        
        if not np.isfinite(correlation_depth) or correlation_depth <= 0:
            raise ValueError("correlation_depth must be finite and positive.")
        if not np.isfinite(timestep) or timestep <= 0:
            raise ValueError("timestep must be finite and positive.")
        self.correlation_depth = float(correlation_depth)
        self.timestep = float(timestep)
        
        # Allow floating-point roundoff at an exact multiple of timestep.
        ratio = self.correlation_depth / self.timestep
        self.files_within_correlation_depth = int(np.floor(np.nextafter(ratio, np.inf)))
        self.structures = deque(maxlen=self.files_within_correlation_depth)
        self.msd_sums = np.zeros(self.files_within_correlation_depth + 1)
        self.msd_counts = np.zeros(self.files_within_correlation_depth + 1, dtype=np.int64)
        self.results = None
        self._nodes = None

    def set_update_structures(self, coords):
        """
        Retain an independent coordinate snapshot in the bounded history.
        """
        self.structures.append(coords.copy())

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
            raise ValueError("MSD requires a nonempty set of finite 3D coordinates.")
        if self._nodes is not None and nodes != self._nodes:
            raise ValueError("MSD atom identities and ordering must remain unchanged.")
        
        self._nodes = nodes
        self.results = None
        
        for lag, previous in enumerate(reversed(self.structures), start=1):
            sample = np.mean(np.sum((coords - previous) ** 2, axis=1))
            self.msd_sums[lag] += sample
            self.msd_counts[lag] += 1
            
        self.set_update_structures(coords)
        return None

    def finalize(self):
        """
        Calculate time-origin averages without writing output.
        """
        lags = np.flatnonzero(self.msd_counts)
        self.results = {
            'lag_time': lags * self.timestep,
            'msd': self.msd_sums[lags] / self.msd_counts[lags],
            'n_origins': self.msd_counts[lags].copy(),
        }
        return self.results

    def render_output(self, data, corr_time):
        self.output_handler.append_row(f"{corr_time:.12g},{data:.12g}")

    def finalize_output(self):
        """Write one averaged result per sampled lag after finalize()."""
        if self.output_handler is None:
            return
        if self.results is None:
            raise RuntimeError("Call finalize() before finalize_output().")
        for lag_time, msd in zip(self.results['lag_time'], self.results['msd']):
            self.render_output(msd, lag_time)
        super().finalize_output(self.header)
