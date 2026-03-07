import numpy as np


class DNAAxisCalculator:

    def __init__(self, args, universe=None, **kwargs):

        if universe is None:
            raise ValueError("DNAAxisCalculator requires an MDAnalysis universe")

        self.args = args
        self.universe = universe
        self.traj = universe.trajectory

        dna = universe.select_atoms("nucleic")

        residues = list(dna.residues)
        residues.sort(key=lambda r: r.resid)

        if len(residues) % 2 != 0:
            raise ValueError("DNA residue count must be even")

        half = len(residues) // 2

        self.strand1 = residues[:half]
        self.strand2 = residues[half:]

        self.axis_per_frame = []

        # smoothing strength
        self.window = 10

    def initialize_output(self):

        self.axis_per_frame = []

    def smooth_axis(self, points):

        points = np.array(points)
        n = len(points)
        w = self.window

        pad = w // 2

        # mirror padding
        padded = np.vstack([
            points[1:pad+1][::-1],
            points,
            points[-pad-1:-1][::-1]
        ])

        smoothed = []

        for i in range(n):
            segment = padded[i:i+w]
            smoothed.append(segment.mean(axis=0))

        return np.array(smoothed)

    def compute_basepair_frame(self, res1, res2):
        """
        Compute base-pair midpoint and normal vector
        """

        p1 = res1.atoms.center_of_mass()
        p2 = res2.atoms.center_of_mass()

        midpoint = 0.5 * (p1 + p2)

        # base pair vector
        x_axis = p2 - p1
        x_axis /= np.linalg.norm(x_axis)

        # approximate base normal using backbone atoms
        a1 = res1.atoms.positions.mean(axis=0)
        a2 = res2.atoms.positions.mean(axis=0)

        y_axis = a2 - a1
        y_axis /= np.linalg.norm(y_axis)

        # base pair normal
        z_axis = np.cross(x_axis, y_axis)
        z_axis /= np.linalg.norm(z_axis)

        return midpoint, z_axis

    def run(self, graph, frame_index):

        ts = self.traj[frame_index]

        midpoints = []
        normals = []

        for i in range(len(self.strand1)):

            res1 = self.strand1[i]
            res2 = self.strand2[-i - 1]

            midpoint, normal = self.compute_basepair_frame(res1, res2)

            midpoints.append(midpoint)
            normals.append(normal)

        midpoints = np.array(midpoints)
        normals = np.array(normals)

        # smooth midpoints
        axis = self.smooth_axis(midpoints)

        # smooth normals
        normals = self.smooth_axis(normals)

        # normalize normals
        normals = normals / np.linalg.norm(normals, axis=1)[:, None]

        # project axis slightly along normals to center helix
        corrected_axis = axis + 0.5 * normals

        self.axis_per_frame.append(corrected_axis)

    def finalize_output(self):

        if len(self.axis_per_frame) == 0:
            return

        filename = "dna_axis.xyz"

        with open(filename, "w") as f:

            for frame in self.axis_per_frame:

                n = len(frame)

                f.write(f"{n}\n")
                f.write("DNA helical axis\n")

                for p in frame:
                    f.write(f"C {p[0]:.5f} {p[1]:.5f} {p[2]:.5f}\n")

        print(f"DNA axis written to {filename}")