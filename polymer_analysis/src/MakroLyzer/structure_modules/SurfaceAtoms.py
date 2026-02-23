import numpy as np
import pytim
import MDAnalysis as mda
import networkx as nx
from collections import defaultdict

from MakroLyzer.structure_modules.structureBase import StructureAnalyzer

class SurfaceAtomsAnalyzer(StructureAnalyzer):
    """
    Analyzer to find the number of nonpolar, polar and ring atoms on the surface of a particle.
    Inherits from StructureAnalyzer.
    """
    
    def __init__(self, output_handler=None, box_size=None, ring_size=6, alpha=5.0):
        """
        Initializes the SurfaceAtomsAnalyzer with the given output handler, box size and alpha parameter.

        Args:
            output_handler (OutputHandler): Handler for writing output.
            box_size (float): The size of the simulation box (assumed cubic).
            ring_size (int|tuple): Ring size filter used for ring identification.
                - int: only this size
                - tuple(min_size, max_size): inclusive range
            alpha (float): The alpha parameter for the gitim algorithm.
        """
        super().__init__(output_handler)
        self.box_size = box_size
        if isinstance(ring_size, int):
            if ring_size < 3:
                raise ValueError("Minimum ring size is 3.")
        elif (
            isinstance(ring_size, tuple)
            and len(ring_size) == 2
            and all(isinstance(v, int) for v in ring_size)
        ):
            min_size, max_size = ring_size
            if min_size < 3:
                raise ValueError("Minimum ring size is 3.")
            if max_size < min_size:
                raise ValueError("Ring-size range must satisfy min <= max.")
        else:
            raise ValueError("ring_size must be an integer or a tuple (min, max).")
        self.ring_size = ring_size
        self.alpha = alpha
        
    def _find_surface_atoms(self, graph):
        """
        We use the gitim algorithm of pytim to identify the surface atoms of the given particle.
    
        Args:
            graph (GraphManager): The graph to modify.
        """
        nodes, coords = graph.get_all_coordinates()
        if len(nodes) == 0:
            return []

        elements = [graph.nodes[n]['element'] for n in nodes]
        coords = np.asarray(coords, dtype=float)

        u = mda.Universe.empty(len(nodes), trajectory=True)
        u.add_TopologyAttr('name', elements)
        u.add_TopologyAttr('type', elements)
        u.add_TopologyAttr('ids', np.arange(1, len(nodes) + 1))
        try:
            u.add_TopologyAttr('elements', elements)
        except Exception:
            pass

        if self.box_size is None:
            min_coords = coords.min(axis=0)
            max_coords = coords.max(axis=0)
            padding = 2.0
            dims = (max_coords - min_coords) + 2.0 * padding
            dims = np.where(dims <= 0, 1.0, dims)
            u.atoms.positions = coords - min_coords + padding
        else:
            if isinstance(self.box_size, (int, float)):
                dims = np.array([self.box_size, self.box_size, self.box_size], dtype=float)
            else:
                dims = np.array(self.box_size, dtype=float)
                if dims.size != 3:
                    dims = np.array([dims[0], dims[0], dims[0]], dtype=float)
            u.atoms.positions = coords

        u.dimensions = [dims[0], dims[1], dims[2], 90.0, 90.0, 90.0]

        interface = pytim.GITIM(
            u,
            group=u.atoms,
            alpha=self.alpha,
            molecular=False,
        )
        if not interface.layers:
            return []

        surface_indices = interface.layers[0].indices
        return [nodes[i] for i in surface_indices]
    
    def _assign_atom_types(self, graph, surface_nodes):
        """
        Assigns atom types (nonpolar, polar, ring) to the surface nodes based on their element type.

        Args:
            graph (GraphManager): The graph containing the nodes.
            surface_nodes (list): List of surface node identifiers.

        Returns:
            dict: A dictionary with counts of each atom type on the surface.
        """
        atom_types = {'nonpolar': 0, 'polar': 0, 'ring': 0}
        
        # find cycles in the graph to identify ring groups
        cycles = nx.cycle_basis(graph)
        cycle_nodes = set()
        for cycle in cycles:
            if isinstance(self.ring_size, int):
                if len(cycle) == self.ring_size:
                    cycle_nodes.update(cycle)
            else:
                min_size, max_size = self.ring_size
                if min_size <= len(cycle) <= max_size:
                    cycle_nodes.update(cycle)
        
        # Assign types
        for node in surface_nodes:
            element = graph.nodes[node]['element']
            if node in cycle_nodes and element == 'C':
                atom_types['ring'] += 1
            elif element in ['C']:
                atom_types['nonpolar'] += 1
            elif element in ['O', 'N', 'S']:
                atom_types['polar'] += 1
                
        return atom_types
    
    def compute(self, graph):
        """
        Compute the number of nonpolar, polar and ring atoms on the surface of the particle.

        Args:
            graph (GraphManager): The graph to analyze.
        Returns:
            dict: A dictionary with counts of each atom type on the surface.
        """
        surface_nodes = self._find_surface_atoms(graph)
        atom_types = self._assign_atom_types(graph, surface_nodes)
        return atom_types
    
    def render_output(self, data, frame_idx):
        """
        Write/Save data for this frame.

        Args:
            data (dict): The computed atom type counts.
            frame_idx (int): Current frame number.
        """
        row = (
            f"{frame_idx},"
            f"{data['nonpolar']},"
            f"{data['polar']},"
            f"{data['ring']}"
        )
        self.output_handler.append_row(row)
        
    def initialize_output(self):
        """
        Initialize output file with header. ('streaming' mode)
        """
        if self.output_handler:
            header = "Frame, Nonpolar Surface Atoms, Polar Surface Atoms, Ring Surface Atoms"
            if self.output_handler.mode == 'streaming':
                self.output_handler.initialize_file(header)
                
    def finalize_output(self):
        """
        Finalize output file (write header and rows - 'collect' mode)
        """
        header = "Frame, Nonpolar Surface Atoms, Polar Surface Atoms, Ring Surface Atoms"
        self.output_handler.finalize(header)
        
    def run(self, graph, frame_idx):
        """
        Execute analysis for a single frame: compute + output

        Args:
            graph : The molecular graph to analyze
            frame_idx (int): Current frame number
        Returns: 
            Analysis Result.
        """        
        result = self.compute(graph)
        if self.output_handler:
            self.render_output(result, frame_idx)
        return result
