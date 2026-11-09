"""Read labeled XYZ templates and explicit bond tables into GraphManager."""

from collections import Counter, defaultdict
from pathlib import Path
import math
import re

from MakroLyzer.graphs import GraphManager


def load_template_graph(
    coordinate_file, bond_file, monomer_type, monomer_id=0, start_index=0,
    *, allowed_elements=None, skip_missing_bonds=False,
):
    """
    Return a fresh graph without inferring bonds from coordinates.

    XYZ rows contain ``label x y z``; labels begin with an element symbol,
    optionally followed by an underscore and site name. Bond tables have a
    header followed by ``label label integer_order`` rows. A repeated label
    expands to all matching atoms when the other endpoint is unique.

    Original labels are stored as ``template_label``; repeated atoms receive
    numbered ``atom_label`` values. Node IDs begin at ``start_index``.
    Missing bond endpoints raise by default; shared tables can explicitly
    enable ``skip_missing_bonds``. ``allowed_elements`` optionally restricts
    element symbols. Coordinates retain their original reference frame.
    """
    if isinstance(start_index, bool) or not isinstance(start_index, int) or start_index < 0:
        raise ValueError("start_index must be a nonnegative integer.")

    path = Path(coordinate_file)
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
        
        if not re.fullmatch(r"[A-Z][a-z]?", element) or (
            allowed_elements is not None and element not in allowed_elements
        ):
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
            monomer_type=monomer_type, x=x, y=y, z=z,
        )
        
        nodes_by_label[label].append(node)

    bond_path = Path(bond_file)
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
        # Shared bond tables may include optional substituents.
        if left not in nodes_by_label or right not in nodes_by_label:
            if skip_missing_bonds:
                continue
            raise ValueError(f"Unknown bond endpoint at row {line_number} in {bond_path}.")
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

