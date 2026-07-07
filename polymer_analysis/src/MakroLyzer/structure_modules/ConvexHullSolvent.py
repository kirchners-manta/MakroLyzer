from collections import deque

import numpy as np
from scipy.spatial import ConvexHull, QhullError

from MakroLyzer import graphs
from MakroLyzer import dictionaries
from MakroLyzer.input_handling import subgraphSelection
from MakroLyzer.structure_modules.structureBase import StructureAnalyzer


class ConvexHullSolventAnalyzer(StructureAnalyzer):
    """
    Analyzer for counting selected solvent molecules inside the convex hull of
    selected particle atoms.
    """

    requires_full_graph = True

    def __init__(self, particle_selection, solvent_selection, box_size=None, output_handler=None, static_topology=False):
        """
        Initialize ConvexHullSolventAnalyzer.

        Args:
            particle_selection (str): Subgraph selection used to build the convex hull.
            solvent_selection (str): Subgraph selection for solvent molecules to count.
            box_size (float | np.ndarray): The size of the periodic box.
            output_handler (OutputHandler): Handler for writing output.
            static_topology (bool): Cache selected nodes after the first frame.
        """
        super().__init__(output_handler)
        self.particle_selection_label = particle_selection
        self.solvent_selection_label = solvent_selection
        self.particle_selection = subgraphSelection.parse_selection_list([particle_selection])
        self.solvent_selection = subgraphSelection.parse_selection_list([solvent_selection])
        self.box_size = box_size
        self.static_topology = static_topology
        self.particle_nodes = None
        self.solvent_nodes = None
     
    @staticmethod
    def _normalize_box_size(box_size):
        box_size = np.asarray(box_size, dtype=float)
        if box_size.ndim == 0:
            box_size = np.repeat(box_size, 3)
        return box_size

    @staticmethod
    def _shift_particle_coordinates(coords: np.ndarray, box_size):
        """
        Shift coordinates to be within the periodic box.

        Args:
            coords (np.ndarray): The coordinates to shift, shape (N, 3) - particle_selection
            box_size (float | np.ndarray): Box length (scalar) or per-axis lengths.

        Returns:
            np.ndarray: A copy of the coordinates with (hopefully) no particles split across the periodic boundary.
        """
        coords   = np.asarray(coords, dtype=float)
        box_size = ConvexHullSolventAnalyzer._normalize_box_size(box_size)
        shifts = graphs.periodic_compacting_shifts(coords, box_size)
        shifted = np.mod(coords + shifts, box_size)

        return shifted, shifts
    
    def _shift_solvent_coordinates(self, coords: np.ndarray, shifts: np.ndarray, box_size):
        """
        Shift solvent coordinates based on the shifts applied to the particle selection.

        Args:
            coords (np.ndarray): The coordinates of the solvent molecules, shape (N, 3).
            shifts (np.ndarray): The shifts applied to the particle selection, shape (3,).
            box_size (float | np.ndarray): Box length (scalar) or per-axis lengths.

        Returns:
            np.ndarray: A copy of the solvent coordinates with shifts applied and wrapped into the periodic box.
        """
        coords   = np.asarray(coords, dtype=float)
        box_size = self._normalize_box_size(box_size)

        shifted = coords.copy()
        # Apply the same shifts to the solvent coordinates and wrap them into the periodic box
        for ax in range(3):
            shifted[:, ax] = np.mod(shifted[:, ax] + shifts[ax], box_size[ax])

        return shifted
    
    def _subgraph_com(self, solvent_subgraph, shifted_coords_by_node):
        """
        Calculate the center of mass for a given solvent subgraph.

        Args:
            solvent_subgraph (GraphManager): The subgraph representing a solvent molecule.
            shifted_coords_by_node (dict): Mapping from node ID to shifted coordinates.
            
        Returns:
            np.ndarray: The center of mass coordinates of the solvent molecule, shape (3,).
        """
        masses = np.array(
            [
                dictionaries.dictMass().get(solvent_subgraph.nodes[node]['element'], 1.0)
                for node in solvent_subgraph.nodes
            ],
            dtype=float,
        )
        # Unwrap the coordinates of the solvent subgraph to account for periodic boundary conditions
        # We need to ensure that the coordinates of the solvent molecule are consistent and not split across periodic boundaries.
        coords_by_node = self._unwrap_subgraph_coordinates(solvent_subgraph, shifted_coords_by_node)
        coords = np.array([coords_by_node[node] for node in solvent_subgraph.nodes], dtype=float)
        total_mass = np.sum(masses)
        if total_mass == 0:
            raise ValueError("Total mass of the solvent molecule is zero, cannot compute center of mass.")
        
        # Calculate the center of mass
        com = np.sum(coords * masses[:, np.newaxis], axis=0) / total_mass
        return np.mod(com, self._normalize_box_size(self.box_size))

    def _unwrap_subgraph_coordinates(self, solvent_subgraph, shifted_coords_by_node):
        nodes = list(solvent_subgraph.nodes)
        if not nodes:
            return {}

        box_size = self._normalize_box_size(self.box_size)

        # Choose one atom as the reference point. Its shifted/wrapped coordinate
        # is kept as-is; all bonded atoms are then placed next to already
        # unwrapped atoms. The chosen root does not matter for the COM, because
        # it only changes the absolute periodic image, not the molecule shape.
        root = nodes[0]
        unwrapped = {root: np.asarray(shifted_coords_by_node[root], dtype=float)}
        queue = deque([root])

        # Breadth-first search follows the molecular connectivity. This is
        # important because a molecule can be split by the periodic boundary
        # even after all atoms were wrapped into the box.
        while queue:
            node = queue.popleft()
            node_wrapped = np.asarray(shifted_coords_by_node[node], dtype=float)
            node_unwrapped = unwrapped[node]

            for neighbor in solvent_subgraph.neighbors(node):
                if neighbor in unwrapped:
                    continue

                # Example: in a box of length 20, bonded atoms at x=19.9 and
                # x=0.1 look 19.8 Å apart in wrapped coordinates, but their
                # minimum-image displacement is only 0.2 Å across the boundary.
                neighbor_wrapped = np.asarray(shifted_coords_by_node[neighbor], dtype=float)
                delta = np.asarray(graphs.min_image_distance(neighbor_wrapped, node_wrapped, box_size))

                # Store a temporary unwrapped coordinate for the neighbor.
                # This does not modify solvent_subgraph or the full graph; it
                # only gives _subgraph_com() a compact molecule for the COM.
                unwrapped[neighbor] = node_unwrapped + delta
                queue.append(neighbor)

        return unwrapped

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
        if self.box_size is None:
            hull_coords = particle_coords
            solvent_coms = [
                solvent_molecule.get_com()
                for solvent_molecule in solvent_graph.get_subgraphs()
            ]
        else:
            hull_coords, shifts = self._shift_particle_coordinates(particle_coords, self.box_size)
            solvent_node_list = list(solvent_nodes)
            solvent_coords = np.asarray([graph.get_coordinates(node) for node in solvent_node_list], dtype=float)
            shifted_solvent_coords = self._shift_solvent_coordinates(solvent_coords, shifts, self.box_size)
            shifted_coords_by_node = dict(zip(solvent_node_list, shifted_solvent_coords))
            solvent_coms = [
                self._subgraph_com(solvent_molecule, shifted_coords_by_node)
                for solvent_molecule in solvent_graph.get_subgraphs()
            ]

        hull, volume = self._compute_hull(hull_coords)

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
