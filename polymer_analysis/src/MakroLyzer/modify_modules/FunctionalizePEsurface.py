import numpy as np
import pytim
import MDAnalysis as mda

from MakroLyzer.modify_modules.structureModifierBase import StructureModifier
from MakroLyzer import dictionaries

class FunctionalizePEsurfaceModifier(StructureModifier):
    """
    A class for functionalizing polyethylene (PE) features. 
    1. Replacing CH2 by CO
    2. Replacing CH2 by NH
    A number between 0 and 100 indicating the percentage of functionalization.
    Inherits fron StructureModifier.
    """
    
    def __init__(self, percentage, func_type, output_handler=None, box_size=None, alpha=5.0):
        """
        Initializes the FunctionalizePE class with the given percentage and functionalization type.

        Args:
            percentage (int): The percentage of functionalization (0-100).
            func_type (str): The type of functionalization ('CO' or 'NH').
        """
        super().__init__(output_handler)
        if func_type not in ['CO']:
            raise ValueError("func_type must be 'CO'")
        if not (0 <= percentage <= 100):
            raise ValueError("percentage must be between 0 and 100")
        
        self.percentage = percentage
        self.func_type = func_type
        self.box_size = box_size
        self.alpha = alpha
        
    def _find_surface_atoms(self, graph):
        """
        We use the gitim algorithm of pytim to identify the surface atoms of the polyethylene particle.
    

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
            #radii_dict=dictionaries.dictVdW(),
        )
        if not interface.layers:
            return []

        surface_indices = interface.layers[0].indices
        return [nodes[i] for i in surface_indices]
        
    def compute(self, graph):
        """
        Functionalizes the polyethylene surface of the given graph.
        """
        surface_nodes = self._find_surface_atoms(graph)
        if not surface_nodes:
            return graph

        surface_graph = graph.subgraph(surface_nodes)
        graph.functionalizePE(self.percentage, self.func_type, surfaceAtoms=surface_graph)
        return graph

    def render_output(self, graph, frame_idx):
        """
        Write graph xyz to file.

        Args:
            graph (GraphManager): The modified graph.
            frame_idx (int): Current frame number.
        """
        NoNodes = graph.number_of_nodes()
        header = f"{NoNodes}\nelement  x         y        z"
        self.output_handler.write_output(header, graph) 
