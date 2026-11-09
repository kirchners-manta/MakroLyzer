from MakroLyzer.builder_modules.BuilderBase import Builder
from MakroLyzer.builder_modules import globalHelpers
from MakroLyzer.graphs import GraphManager
from pathlib import Path
from collections import Counter, defaultdict
import math

class LigninBuilder(Builder):
    """
    Builds a lignin strand corresponding to the users given sequence of monomers and likages.
    Inherits from Builder.
    """
    
    def __init__(self, recipe_file, output_handler=None):
        super().__init__(output_handler)
        self.recipe = globalHelpers.read_recipe(recipe_file) # globally for all builders
        
    SUPPORTED_MONOMERS = frozenset({'H', 'G', 'S'})
    SUPPORTED_LINKERS = frozenset({
        'beta-O-4', 'beta-5', 'beta-beta', '5-5', '4-O-beta',
        '5-O-4', 'beta-1',
    })

    def load_monomer(self, monomer, monomer_id=0, start_index=0):
        """Load a fresh GraphManager with explicit template connectivity.

        Node IDs start at ``start_index``. When merging monomers, pass an
        unused index range and a distinct ``monomer_id``. ``template_label``
        preserves source labels; ``atom_label`` distinguishes repeated atoms.
        Coordinates are in the template's original reference frame.
        """
        if not isinstance(monomer, str) or monomer not in self.SUPPORTED_MONOMERS:
            raise ValueError(f"Unsupported monomer: {monomer!r}.")
        if isinstance(start_index, bool) or not isinstance(start_index, int) or start_index < 0:
            raise ValueError("start_index must be a nonnegative integer.")

        directory = Path(__file__).resolve().parent / 'buildingBlocks' / 'lignin'
        path = directory / f'{monomer}-topo.xyz'
        lines = path.read_text(encoding='utf-8').splitlines()
        if len(lines) < 2:
            raise ValueError(f"Incomplete XYZ template: {path}.")
        try:
            atom_count = int(lines[0])
        except ValueError as exc:
            raise ValueError(f"Invalid atom count in {path}.") from exc
        rows = [line.split() for line in lines[2:] if line.strip()]
        if atom_count < 1 or len(rows) != atom_count:
            raise ValueError(f"Atom count does not match template rows in {path}.")

        counts = Counter(row[0] for row in rows)
        occurrences = Counter()
        nodes_by_label = defaultdict(list)
        graph = GraphManager()
        for offset, row in enumerate(rows):
            if len(row) != 4:
                raise ValueError(f"Invalid atom row {offset + 3} in {path}.")
            label = row[0]
            element = label.split('_', 1)[0]
            if element not in {'C', 'H', 'O'}:
                raise ValueError(f"Unsupported atom label {label!r} in {path}.")
            try:
                x, y, z = map(float, row[1:])
            except ValueError as exc:
                raise ValueError(f"Invalid coordinates at row {offset + 3} in {path}.") from exc
            if not all(math.isfinite(value) for value in (x, y, z)):
                raise ValueError(f"Nonfinite coordinates at row {offset + 3} in {path}.")
            occurrences[label] += 1
            atom_label = f'{label}_{occurrences[label]}' if counts[label] > 1 else label
            node = start_index + offset
            graph.add_node(
                node, index=node, element=element, atom_label=atom_label,
                template_label=label, monomer_id=monomer_id,
                monomer_type=monomer, x=x, y=y, z=z,
            )
            nodes_by_label[label].append(node)

        bond_path = directory / 'bonding-info.txt'
        for line_number, line in enumerate(bond_path.read_text(encoding='utf-8').splitlines()[1:], 2):
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(f"Invalid bond row {line_number} in {bond_path}.")
            left, right, order = fields
            order = int(order)
            if order not in {1, 2, 3}:
                raise ValueError(f"Invalid bond order at row {line_number} in {bond_path}.")
            # The shared table includes substituents absent from some monomers.
            if left not in nodes_by_label or right not in nodes_by_label:
                continue
            if len(nodes_by_label[left]) > 1 and len(nodes_by_label[right]) > 1:
                raise ValueError(f"Ambiguous bond at row {line_number} in {bond_path}.")
            for source in nodes_by_label[left]:
                for target in nodes_by_label[right]:
                    graph.add_edge(source, target, bond_order=order)

        for node, degree in graph.degree():
            if degree == 0:
                raise ValueError(f"Unbonded atom {graph.nodes[node]['atom_label']} in {path}.")
            graph.nodes[node]['degree'] = degree
        return graph

    def _connect_monomers(self, chain, linker, monomer):
        raise NotImplementedError("Monomer connection is not implemented yet.")
        
    @staticmethod
    def _forbidden_combinations():
        """Return ordered (left monomer, linker, right monomer) exclusions.

        Map each forbidden triple to a human-readable reason. Rules are
        directional; add the reverse triple explicitly when needed.
        This is a scaffold for confirmed linkage rules, not a complete
        chemical compatibility check.
        """
        return {
            # ('left monomer', 'linker', 'right monomer'): 'Reason',
        }
        
    def validate_recipe(self):
        """Validate a linear recipe, raising ValueError for invalid input."""
        if not isinstance(self.recipe, dict):
            raise ValueError("Recipe must be a dictionary.")

        for key in ('monomer', 'linker'):
            if key not in self.recipe:
                raise ValueError(f"Recipe is missing required field '{key}'.")
            if not isinstance(self.recipe[key], list):
                raise ValueError(f"Recipe field '{key}' must be a list.")

        monomers = self.recipe['monomer']
        linkers = self.recipe['linker']
        if not monomers:
            raise ValueError("Recipe must contain at least one monomer.")
        if len(monomers) != len(linkers) + 1:
            raise ValueError(
                "A linear recipe requires exactly one fewer linker than monomers "
                f"(got {len(monomers)} monomers and {len(linkers)} linkers)."
            )

        for key, values, supported in (
            ('monomer', monomers, self.SUPPORTED_MONOMERS),
            ('linker', linkers, self.SUPPORTED_LINKERS),
        ):
            for i, value in enumerate(values):
                if not isinstance(value, str) or value not in supported:
                    raise ValueError(
                        f"Unsupported {key} at index {i}: {value!r}. "
                        f"Supported values: {', '.join(sorted(supported))}."
                    )

        forbidden = self._forbidden_combinations()
        for i, linker in enumerate(linkers):
            combination = (monomers[i], linker, monomers[i + 1])
            if combination in forbidden:
                raise ValueError(
                    f"Forbidden combination at linker index {i}: "
                    f"{combination!r}. {forbidden[combination]}"
                )
        
    def build_structure(self):
        """
        1) Load first monomer and call it chain
        2) For each -linker-monomer tuple:
            - load monomer
            - connect monober to chain via the linker
            - energy optimization        
        """
        chain = self.load_monomer(self.recipe["monomer"][0])
        for i, linker in enumerate(self.recipe["linker"]):
            monomer = self.load_monomer(
                self.recipe["monomer"][i + 1], monomer_id=i + 1,
                start_index=max(chain.nodes) + 1,
            )
            chain = self._connect_monomers(chain, linker, monomer)
            chain = self.minimize_energy(chain)
            
        return chain
    
    def render_output(self, data):
        """Write the built lignin graph as XYZ and return the output path."""
        if self.output_handler is None:
            return None
        header = f"{data.number_of_nodes()}\nLignin structure"
        return self.output_handler.write_output(header, data)
    
