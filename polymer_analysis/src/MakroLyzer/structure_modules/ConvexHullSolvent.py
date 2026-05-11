import numpy as np
from scipy.spatial import ConvexHull, QhullError

from MakroLyzer import graphs
from MakroLyzer.input_handling import subgraphSelection
from MakroLyzer.structure_modules.structureBase import StructureAnalyzer


class ConvexHullSolventAnalyzer(StructureAnalyzer):
    """
    Analyzer for counting selected solvent molecules inside the convex hull of
    selected particle atoms.
    """

    requires_full_graph = True

    def __init__(self, particle_selection, solvent_selection, output_handler=None, static_topology=False):
        """
        Initialize ConvexHullSolventAnalyzer.

        Args:
            particle_selection (str): Subgraph selection used to build the convex hull.
            solvent_selection (str): Subgraph selection for solvent molecules to count.
            output_handler (OutputHandler): Handler for writing output.
            static_topology (bool): Cache selected nodes after the first frame.
        """
        super().__init__(output_handler)
        self.particle_selection_label = particle_selection
        self.solvent_selection_label = solvent_selection
        self.particle_selection = subgraphSelection.parse_selection_list([particle_selection])
        self.solvent_selection = subgraphSelection.parse_selection_list([solvent_selection])
        self.static_topology = static_topology
        self.particle_nodes = None
        self.solvent_nodes = None

    def initialize_output(self):
        """
        Initialize output file with header. ('streaming' mode)
        """
        if self.output_handler:
            header = (
                "Frame, Particle Selection, Solvent Selection, Convex Hull - Volume / Å³, "
                "Number of Solvent Molecules Inside"
            )
            if self.output_handler.mode == 'streaming':
                self.output_handler.initialize_file(header)

    def _get_selection_nodes(self, graph):
        if self.particle_nodes is None or self.solvent_nodes is None or not self.static_topology:
            self.particle_nodes = graph.select_subgraph_nodes(self.particle_selection)
            self.solvent_nodes = graph.select_subgraph_nodes(self.solvent_selection)
        return self.particle_nodes, self.solvent_nodes

    @staticmethod
    def _compute_hull(coords):
        if len(coords) < 4:
            return None, 0.0

        try:
            hull = ConvexHull(coords)
        except QhullError:
            return None, 0.0
        return hull, float(hull.volume)

    @staticmethod
    def _points_inside_hull(hull, points, tolerance=1e-12):
        if hull is None or len(points) == 0:
            return np.zeros(len(points), dtype=bool)

        # We obtain the normal vectors and offsets of the hull facets.
        # All normal vectors point outward from the hull, so for a point to be inside the hull,
        # the dot product of the point with the normal vector plus the offset should be <= 0 for all facets.
        equations = hull.equations
        # equations[:, :-1].T : Normal vectors of the hull facets
        # equations[:, -1] : Offsets of the hull facets
        return np.all(np.dot(points, equations[:, :-1].T) + equations[:, -1] <= tolerance, axis=1)

    def compute(self, graph):
        """
        Calculate the selected particle's convex hull and count selected solvent
        molecules whose center of mass lies inside the hull.

        Args:
            graph (GraphManager): Full graph to analyze.

        Returns:
            tuple: particle selection, solvent selection, hull volume, solvent count.
        """
        particle_nodes, solvent_nodes = self._get_selection_nodes(graph)
        solvent_graph = graphs.GraphManager(graph.subgraph(solvent_nodes))

        particle_coords = np.asarray([graph.get_coordinates(node) for node in particle_nodes], dtype=float)
        hull, volume = self._compute_hull(particle_coords)

        solvent_coms = []
        for solvent_molecule in solvent_graph.get_subgraphs():
            solvent_coms.append(solvent_molecule.get_com())

        inside = self._points_inside_hull(hull, np.asarray(solvent_coms, dtype=float))
        count = int(np.count_nonzero(inside))

        return (
            self.particle_selection_label,
            self.solvent_selection_label,
            volume,
            count,
        )

    def render_output(self, data, frame_idx):
        """
        Write/Save data for this frame.

        Args:
            data (tuple): The computed convex hull solvent data.
            frame_idx (int): Current frame number.
        """
        particle_selection, solvent_selection, volume, count = data
        row = (
            f"{frame_idx},"
            f"{particle_selection},"
            f"{solvent_selection},"
            f"{volume:.3f},"
            f"{count}"
        )
        self.output_handler.append_row(row)

    def finalize_output(self):
        """
        Finalize output file (write header and rows - 'collect' mode).
        """
        header = (
            "Frame, Particle Selection, Solvent Selection, Convex Hull - Volume / Å³, "
            "Number of Solvent Molecules Inside"
        )
        super().finalize_output(header)
