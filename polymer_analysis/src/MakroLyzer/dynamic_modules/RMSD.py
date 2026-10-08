import numpy as np

from MakroLyzer.dynamic_modules.dynamicBase import DynamicAnalyzer


class RMSDAnalyzer(DynamicAnalyzer):
    """Atom-averaged RMSD relative to the first analyzed frame.

    Return raw RMSD and RMSD after subtracting each frame's mass-weighted
    center of mass. Both RMSDs weight atoms equally; no rotational alignment
    or averaging over time origins is performed. Coordinates must use
    consistent periodic images (the dynamic driver unwraps the trajectory).
    """

    header = "Frame, RMSD / Å, COM-corrected RMSD / Å"

    def __init__(self, output_handler=None):
        super().__init__(output_handler)
        self.refstruc = None
        self.refstruc_com = None
        self._nodes = None

    def set_first_structure(self, coords):
        """
        Copy reference coordinates to allow in-place graph updates.
        """
        self.refstruc = np.asarray(coords, dtype=float).copy()

    def get_refstruc_com(self, graph):
        """
        Copy the reference center of mass.
        """
        self.refstruc_com = np.asarray(graph.get_com(), dtype=float).copy()

    def initialize_output(self):
        if self.output_handler is not None and self.output_handler.mode == 'streaming':
            self.output_handler.initialize_file(self.header)

    def compute(self, graph):
        """
        Return (raw RMSD, COM-corrected RMSD) for the current frame.
        """
        nodes, coords = graph.get_all_coordinates()
        nodes = tuple(nodes)
        coords = np.asarray(coords, dtype=float)
        if not nodes or coords.shape != (len(nodes), 3) or not np.all(np.isfinite(coords)):
            raise ValueError("RMSD requires a nonempty set of finite 3D coordinates.")
        if self._nodes is not None and nodes != self._nodes:
            raise ValueError("RMSD atom identities and ordering must remain unchanged.")

        com = np.asarray(graph.get_com(), dtype=float)
        if com.shape != (3,) or not np.all(np.isfinite(com)):
            raise ValueError("RMSD requires a finite 3D center of mass.")
        if self.refstruc is None:
            self.set_first_structure(coords)
        if self.refstruc_com is None:
            self.refstruc_com = com.copy()
        self._nodes = nodes

        displacement = coords - self.refstruc
        corrected_displacement = (coords - com) - (self.refstruc - self.refstruc_com)
        rmsd = np.sqrt(np.mean(np.sum(displacement ** 2, axis=1)))
        corrected_rmsd = np.sqrt(np.mean(np.sum(corrected_displacement ** 2, axis=1)))
        return rmsd, corrected_rmsd

    def render_output(self, data, frame_idx):
        """
        Write both RMSDs for this frame, in angstroms.
        """
        rmsd, corrected_rmsd = data
        self.output_handler.append_row(f"{frame_idx},{rmsd:.3f},{corrected_rmsd:.3f}")

    def finalize_output(self):
        super().finalize_output(self.header)
