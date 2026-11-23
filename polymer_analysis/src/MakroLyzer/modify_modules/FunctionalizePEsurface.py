import numpy as np
import pytim
import MDAnalysis as mda

from MakroLyzer.modify_modules.FunctionalizePE import write_functionalization_output
from MakroLyzer.modify_modules.structureModifierBase import StructureModifier

class FunctionalizePEsurfaceModifier(StructureModifier):
    """Randomly functionalize carbons in the outermost GITIM layer of PE.

    Targets and report percentages use the original surface carbon count.
    Surface chain ends contribute to that count and are excluded from placement.
    """

    def __init__(self, percentage, func_type, output_handler=None, box_size=None, alpha=5.0, neighbor_exclusion=1,
                 *, requests=None, original_input=None):
        """Configure functionalization and optional output.

        Args:
            percentage (float): Percentage of original surface carbons to
                functionalize (0-100).
            func_type (str): One of CO, COH, ester, or amide.
            output_handler (ModifyOutputHandler, optional): Writer for the XYZ output.
            box_size (float or sequence, optional): Orthorhombic box lengths in
                angstroms. A scalar sets all three lengths. If omitted, the box
                is estimated from the coordinate bounds.
            alpha (float): GITIM probe radius in angstroms.
            neighbor_exclusion (int): Inclusive exclusion distance in C-C
                bonds along the full backbone. Defaults to one.
            requests (list[tuple[str, float]], optional): Ordered
                (group type, percentage) pairs for mixed random placement.
                Each group reserves sites in input order using the same
                exclusion distance. When supplied, these pairs set all group
                targets. The single-group arguments still undergo validation.
            original_input (str, optional): Raw user input to record in the
                companion .info.txt report.
        """
        super().__init__(output_handler)
        if func_type not in ['CO', 'COH', 'ester', 'amide']:
            raise ValueError("func_type must be 'CO', 'COH','ester' or 'amide'")
        if not (0 <= percentage <= 100):
            raise ValueError("percentage must be between 0 and 100")
        
        self.func_type = func_type
        self.percentage = percentage
        self.box_size = box_size
        self.alpha = alpha
        self.neighbor_exclusion = neighbor_exclusion
        self.requests = requests
        self.original_input = original_input
        self.report = None
        
    def _find_surface_atoms(self, graph):
        """Return graph node IDs from the outermost GITIM layer.

        Pass all atoms to PyTIM and map the detected layer indices back to
        the original graph. Return an empty list when there is no layer.
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
        
    def compute(self, graph):
        """Functionalize the detected surface carbons in place.

        Store the counts and percentages in `self.report` and return the
        input graph with the modifications applied.
        """
        surface_nodes = self._find_surface_atoms(graph)
        # Pass an empty subgraph when GITIM detects no surface atoms so the report
        # records zero surface carbons.
        surface_graph = graph.subgraph(surface_nodes)
        options = {} if self.requests is None else {'requests': self.requests}
        self.report = graph.functionalizePE(
            mode="random",
            func_type=self.func_type,
            percentage=self.percentage,
            neighbor_exclusion=self.neighbor_exclusion,
            surface_atoms=surface_graph,
            **options,
        )
        return graph

    def render_output(self, graph, frame_idx):
        """Write the modified XYZ and its accompanying .info.txt report.

        Args:
            graph (GraphManager): Graph returned by compute().
            frame_idx (int): Frame index supplied by StructureModifier.run().
                The output writer chooses the filename.
        """
        write_functionalization_output(self.output_handler, graph, self.report, self.original_input)
