"""Lignin linkage construction and geometry hooks, independent of the builder."""

import numpy as np
from scipy.spatial.transform import Rotation
from MakroLyzer import dictionaries


def _find_attachment_atom(graph, template_label, monomer_id):
    """Find the atom of the monomer in the chain to connect the new monomer to"""
    # Labels recur in different monomers, so both label and monomer ID are needed.
    matches = [
        node for node, data in graph.nodes(data=True)
        if data.get('template_label') == template_label
        and data.get('monomer_id') == monomer_id
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one {template_label} in monomer {monomer_id}; found {len(matches)}.")
    return matches[0]

# -----------------------------------------------------------------------------------------------
# Beta-O-4 --------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------

def connect_beta_O_4(chain, monomer, chain_monomer_id):
    """Scaffold: chain beta carbon connects to incoming monomer O-4.

    Find corresponding two atoms, place new monomer at the right position of the chain and 
    add/remove edges.
    """    
    
    if set(chain).intersection(monomer):
        raise ValueError("Chain and incoming monomer must have disjoint node IDs.")
    chain = chain.copy()
    monomer = monomer.copy()
    # The incoming graph represents one new monomer, with an identity not
    # already present in the chain. Node IDs and monomer IDs serve different roles.
    incoming_ids = {data['monomer_id'] for _, data in monomer.nodes(data=True)}
    if len(incoming_ids) != 1:
        raise ValueError("Incoming graph must contain exactly one monomer.")
    incoming_id = incoming_ids.pop()  # Extract the single ID from the set.

    # beta belongs to the preceding recipe monomer.
    # O_4 belongs to the incoming monomer. 
    # These variables hold node IDs, not atom data.
    beta_atom = _find_attachment_atom(chain, 'C_b', chain_monomer_id)
    oxygen_atom = _find_attachment_atom(monomer, 'O_4', incoming_id)
    
    # Place new monomer at the correct position near the existing chain
    place_beta_O_4(chain, monomer, beta_atom, oxygen_atom)
    edit_beta_O_4_topology(chain, monomer, beta_atom, oxygen_atom)
    
    # Add nodes and edges from the incomming monomer to the chain
    result = chain.copy()
    result.add_nodes_from(monomer.nodes(data=True))
    result.add_edges_from(monomer.edges(data=True))
    # Add the new beta-O-4 bond connecting chain and new monomer
    result.add_edge(beta_atom, oxygen_atom, bond_order=1, linkage='beta-O-4')
    
    # update degree
    for node, degree in result.degree():
        result.nodes[node]['degree'] = degree
        
    return result

def place_beta_O_4(chain, monomer, beta_atom, oxygen_atom):
    """
    Place incoming monomer with its O_4 above C_b along the six-ring 
    normal of the "C_b" monomer.

    Rigidly rotate O_4 - C_4 onto that normal, giving C_b - O_4 - C_4
    an initial angle of 180 degrees. 
    Ring label order defines the positive normal
    -> This initial placement does not guarantee absence of clashes.
    """
    chain_id = chain.nodes[beta_atom]['monomer_id']
    incoming_id = monomer.nodes[oxygen_atom]['monomer_id']
    
    # Get ring-normal
    ring = np.array([
        _coordinates(chain, _find_attachment_atom(chain, f'C_{i}', chain_id))
        for i in range(1, 7)
    ])
    ## singular values of a matrix are the square roots of the eigenvalues of A^T*A
    ## tell how much the points spread out along each principle axis in V
    ##          - Largest: maximum spread (along length of the plane)
    ##          - Smallest: minimum spread (thickness of plane)
    _, singular_values, axes = np.linalg.svd(ring - ring.mean(axis=0))
    if singular_values[1] < 1e-10:
        raise ValueError("Cannot define a plane from a degenerate six-ring.")
    normal = axes[-1] # smallest
    reference = np.cross(ring[1] - ring[0], ring[2] - ring[0])
    if np.dot(normal, reference) < 0:
        normal = -normal

    # get C-O bond vector of incomming monomer
    carbon = _find_attachment_atom(monomer, 'C_4', incoming_id)
    origin = _coordinates(monomer, oxygen_atom)
    direction = _coordinates(monomer, carbon) - origin
    length = np.linalg.norm(direction)
    direction /= length
    
    # get vecor that stands perpendiculat on ring-normal and C-O bond vector
    # and rotate C-O bond vector around that 
    axis = np.cross(direction, normal)
    sine = np.linalg.norm(axis)
    cosine = np.clip(np.dot(direction, normal), -1.0, 1.0)
    if sine > 1e-10:
        rotation = Rotation.from_rotvec(axis / sine * np.arctan2(sine, cosine))
    elif cosine < 0:
        # Antiparallel vectors need a half-turn around any perpendicular axis.
        axis = np.cross(direction, np.eye(3)[np.argmin(abs(direction))])
        rotation = Rotation.from_rotvec(axis / np.linalg.norm(axis) * np.pi)
    else:
        rotation = Rotation.identity()

    # update coords of incomming monomer
    bond_length = dictionaries.dictBondDistances()[('C-O', 'ether')]
    target = _coordinates(chain, beta_atom) + bond_length * normal
    nodes = list(monomer)
    positions = np.array([_coordinates(monomer, node) for node in nodes])
    positions = rotation.apply(positions - origin) + target
    for node, position in zip(nodes, positions):
        monomer.nodes[node].update(zip(('x', 'y', 'z'), map(float, position)))


def _coordinates(graph, node):
    """Read a finite coordinate vector from a graph node."""
    position = np.array([graph.nodes[node][key] for key in ('x', 'y', 'z')], dtype=float)
    if not np.isfinite(position).all():
        raise ValueError(f"Nonfinite coordinates for atom {node}.")
    return position

def edit_beta_O_4_topology(chain, monomer, beta_atom, oxygen_atom):
    """Apply the specified beta-side and incoming O_4-side edits in place.

    Reduce C_a - C_b to order 1 and replace H_a with O_a - H_Oa.
    Add the new hydroxyl along the former C_a - H_a direction.
    Attachment atoms are selected by connect_beta_O_4; template bonds are trusted.
    """
    # Restrict alpha-side edits to the same monomer as the selected beta carbon.
    monomer_id = chain.nodes[beta_atom]['monomer_id']
    alpha_atom = _find_attachment_atom(chain, 'C_a', monomer_id)
    alpha_hydrogen = _find_attachment_atom(chain, 'H_a', monomer_id)
    
    # Resolve the incoming phenolic hydrogen within the selected O_4 monomer.
    incoming_id = monomer.nodes[oxygen_atom]['monomer_id']
    phenolic_hydrogen = _find_attachment_atom(monomer, 'H_O4', incoming_id)
    
    # Require the expected starting bond before changing the topology.
    if (
        not chain.has_edge(alpha_atom, beta_atom)
        or chain.edges[alpha_atom, beta_atom].get('bond_order') != 2
    ):
        raise ValueError("Expected a C_a--C_b double bond before beta-O-4 edits.")

    # Get coordinates for OH, which will be added to alpha carbon of "old" monomer
    oxygen_id = max(max(chain.nodes), max(monomer.nodes)) + 1
    hydrogen_id = oxygen_id + 1
    monomer_type = chain.nodes[alpha_atom]['monomer_type']
    alpha_position = _coordinates(chain, alpha_atom)
    hydroxyl_direction = _coordinates(chain, alpha_hydrogen) - alpha_position
    direction_length = np.linalg.norm(hydroxyl_direction)
    if direction_length < 1e-10:
        raise ValueError("C_a and H_a coordinates coincide.")
    hydroxyl_direction /= direction_length
    bond_distances = dictionaries.dictBondDistances()
    carbon_oxygen_distance = bond_distances[('C-O', 'alcohol')]
    oxygen_hydrogen_distance = bond_distances[('O-H', 'alcohol')]
    oxygen_position = alpha_position + carbon_oxygen_distance * hydroxyl_direction
    hydrogen_position = oxygen_position + oxygen_hydrogen_distance * hydroxyl_direction
    
    # Adjust bond order between alpha and beta carbons of "old" monomer
    chain.edges[alpha_atom, beta_atom]['bond_order'] = 1
    
    # Remove the alpha hydrogen from the "old" monomer
    chain.remove_node(alpha_hydrogen)
    
    # Remove the H_O4 atom from the incomming monomer
    monomer.remove_node(phenolic_hydrogen)
    
    # Add OH nodes and edges to "old" monomer
    for node, label, element, position in (
        (oxygen_id, 'O_a', 'O', oxygen_position),
        (hydrogen_id, 'H_Oa', 'H', hydrogen_position),
    ):
        chain.add_node(
            node, index=node, element=element, atom_label=label,
            template_label=label, monomer_id=monomer_id,
            monomer_type=monomer_type,
            x=float(position[0]), y=float(position[1]), z=float(position[2]),
        )
    # Form the alpha hydroxyl group: C_a - O_a - H_Oa.
    chain.add_edge(alpha_atom, oxygen_id, bond_order=1)
    chain.add_edge(oxygen_id, hydrogen_id, bond_order=1)
    for node, degree in chain.degree():
        chain.nodes[node]['degree'] = degree
    for node, degree in monomer.degree():
        monomer.nodes[node]['degree'] = degree

