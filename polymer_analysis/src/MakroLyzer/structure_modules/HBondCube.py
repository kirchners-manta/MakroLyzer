import numpy as np

from MakroLyzer.structure_modules.structureBase import StructureAnalyzer


ANGSTROM_TO_BOHR = 1.8897259886


class HBondCubeAnalyzer(StructureAnalyzer):
    """
    Analyzer that records H-bond counts on a 3D grid and writes a .cube file.
    """

    def __init__(self, cutoffs, box_size, no_cells_per_dim, output_handler=None):
        super().__init__(output_handler)
        self.cutoffs = cutoffs
        self.box_size = np.asarray(box_size, dtype=float)
        self.no_cells_per_dim = np.asarray(no_cells_per_dim, dtype=int)
        self.cell_size = self.box_size / self.no_cells_per_dim
        self.grid = np.zeros(tuple(self.no_cells_per_dim), dtype=np.int64)
        self._n_frames = 0

    def initialize_output(self):
        # No streaming output for cube fields.
        pass

    def _grid_index(self, position):
        wrapped = np.mod(np.asarray(position, dtype=float), self.box_size)
        idx = np.floor(wrapped / self.cell_size).astype(int)
        idx = np.clip(idx, 0, self.no_cells_per_dim - 1)
        return tuple(idx.tolist())

    def compute(self, graph):
        frame_hits = 0
        for element_type, h_acceptor_dist, donor_acceptor_dist, angle_cut in self.cutoffs:
            hbonds = graph.get_hbonds(element_type, h_acceptor_dist, donor_acceptor_dist, angle_cut)
            for node_h, _ in hbonds:
                idx = self._grid_index(graph.get_coordinates(node_h))
                self.grid[idx] += 1
                frame_hits += 1
        self._n_frames += 1
        return frame_hits

    def render_output(self, data, frame_idx):
        # Intentionally no per-frame text output.
        return None

    def _cube_text(self):
        nx, ny, nz = self.no_cells_per_dim
        lx, ly, lz = self.box_size
        origin = np.array([0.0, 0.0, 0.0]) * ANGSTROM_TO_BOHR
        vx = np.array([lx / nx, 0.0, 0.0]) * ANGSTROM_TO_BOHR # voxel vector
        vy = np.array([0.0, ly / ny, 0.0]) * ANGSTROM_TO_BOHR # voxel vector
        vz = np.array([0.0, 0.0, lz / nz]) * ANGSTROM_TO_BOHR # voxel vector

        lines = [
            "MakroLyzer HBond cube",
            "H-bond counts deposited at hydrogen positions",
            f"  1 {origin[0]: .6f} {origin[1]: .6f} {origin[2]: .6f}",
            f"{nx:4d} {vx[0]: .6f} {vx[1]: .6f} {vx[2]: .6f}",
            f"{ny:4d} {vy[0]: .6f} {vy[1]: .6f} {vy[2]: .6f}",
            f"{nz:4d} {vz[0]: .6f} {vz[1]: .6f} {vz[2]: .6f}",
            "  0  0.000000  0.000000  0.000000  0.000000",
        ]

        flat = self.grid.ravel(order="C")
        for start in range(0, flat.size, 6):
            chunk = flat[start:start + 6]
            lines.append(" ".join(f"{val: .5e}" for val in chunk))
        return "\n".join(lines) + "\n"

    def finalize_output(self):
        if not self.output_handler:
            return
        self.output_handler.write_raw(self._cube_text())
