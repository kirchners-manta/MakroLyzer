"""Lignin linkage construction and geometry hooks, independent of the builder."""

import numpy as np
from scipy.spatial.transform import Rotation
from scipy.spatial import cKDTree
import networkx as nx
from MakroLyzer import dictionaries


# -----------------------------------------------------------------------------------------------
# Generell --------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------

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

def clash_score(fixed_graph, moving_graph, excluded_pairs=(), scale=0.75):
    """Sum squared cross-graph overlaps using a KD-tree (angstroms squared).

    A pair clashes below scale * (vdW_radius_a + vdW_radius_b). Exclusions
    are ordered (fixed_node_id, moving_node_id) pairs supplied by the linker;
    excluding one pair does not exempt either atom from other contacts.
    Only cross-graph contacts are scored. Neither input is modified.
    """
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Clash scale must be finite and positive.")
    if not fixed_graph or not moving_graph:
        return 0.0
    excluded_pairs = set(excluded_pairs)
    radii = dictionaries.dictVdW()
    # fixed_graph is the existing chain; moving_graph is the candidate incoming
    # monomer. These names describe their roles, not left/right in physical space.
    # Keep node IDs in lists: array indices i/j are not necessarily graph IDs.
    fixed_nodes = list(fixed_graph)
    moving_nodes = list(moving_graph)
    fixed_positions = np.array([_coordinates(fixed_graph, n) for n in fixed_nodes])
    moving_positions = np.array([_coordinates(moving_graph, n) for n in moving_nodes])
    fixed_radii = np.array([radii[fixed_graph.nodes[n]['element']] for n in fixed_nodes])
    moving_radii = np.array([radii[moving_graph.nodes[n]['element']] for n in moving_nodes])
    # Index chain coordinates spatially so distant atoms need no pairwise check.
    tree = cKDTree(fixed_positions)
    # A conservative per-atom cutoff finds every possible overlap. Actual
    # element-specific thresholds are applied only to the nearby candidates.
    neighbors = tree.query_ball_point(
        moving_positions, scale * (moving_radii + fixed_radii.max()),
    )
    score = 0.0
    # neighbors[j] contains ARRAY INDICES of chain atoms near incoming atom j.
    # Convert those indices to node IDs when checking linker-specific exemptions.
    for j, nearby in enumerate(neighbors):
        indices = np.array([
            i for i in nearby
            if (fixed_nodes[i], moving_nodes[j]) not in excluded_pairs
        ], dtype=int)
        if not len(indices):
            continue
        # Compare each nearby pair with its own element-dependent threshold.
        # Positive overlap = too close; negative values are clipped to zero.
        distances = np.linalg.norm(fixed_positions[indices] - moving_positions[j], axis=1)
        overlaps = np.maximum(scale * (fixed_radii[indices] + moving_radii[j]) - distances, 0)
        score += float(np.dot(overlaps, overlaps))  # Sum of squared overlaps, not energy.
    return score


def _coordinates(graph, node):
    """Read a finite coordinate vector from a graph node."""
    position = np.array([graph.nodes[node][key] for key in ('x', 'y', 'z')], dtype=float)
    if not np.isfinite(position).all():
        raise ValueError(f"Nonfinite coordinates for atom {node}.")
    return position

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
    Try both ring-normal signs and 20-degree rotations about the new bond.
    Choose the least-overlapping candidate, including the new alpha hydroxyl.
    This is an initial placement search, not a guarantee of clash-free geometry.
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
    normal = axes[-1] # Axis with the smallest spread: perpendicular to the ring.
    # SVD may return either sign. Ring atom order makes the initial sign repeatable;
    # the placement search will still consider both sides of the plane.
    reference = np.cross(ring[1] - ring[0], ring[2] - ring[0])
    if np.dot(normal, reference) < 0:
        normal = -normal

    # get C-O bond vector of incomming monomer
    carbon = _find_attachment_atom(monomer, 'C_4', incoming_id)
    origin = _coordinates(monomer, oxygen_atom)
    direction = _coordinates(monomer, carbon) - origin
    length = np.linalg.norm(direction)
    
    direction /= length
    bond_length = dictionaries.dictBondDistances()[('C-O', 'ether')]
    nodes = list(monomer)
    # Shift coordinates to O_4 as origin so rotations pivot about that atom.
    # Preserve this starting array so trials do not accumulate rotation errors.
    original = np.array([_coordinates(monomer, node) for node in nodes]) - origin
    best_score = float('inf')
    best_positions = None

    # Outer loop: attach above or below the preceding monomer's ring plane, 
    #             above the C_b of the "old" monomer
    for sign in (1, -1):
        outward = sign * normal
        # Rotate O_4->C_4 onto the chosen outward direction. The cross product
        # gives the rotation axis; atan2(sine, cosine) gives its angle in radians.
        axis = np.cross(direction, outward)
        sine = np.linalg.norm(axis)
        cosine = np.clip(np.dot(direction, outward), -1.0, 1.0)
        if sine > 1e-10:
            rotation = Rotation.from_rotvec(axis / sine * np.arctan2(sine, cosine))
        elif cosine < 0:
            axis = np.cross(direction, np.eye(3)[np.argmin(abs(direction))])
            rotation = Rotation.from_rotvec(axis / np.linalg.norm(axis) * np.pi)
        else:
            rotation = Rotation.identity()  # Already aligned.
        aligned = rotation.apply(original)
        # Translation puts O_4 one dictionary bond length away from C_b.
        target = _coordinates(chain, beta_atom) + bond_length * outward
        
        # Inner loop: twist the aligned monomer about the new C_b--O_4 axis. (Like a spindle)
        # This changes its orientation without changing the attachment distance.
        for angle in range(0, 360, 20):
            twist = Rotation.from_rotvec(outward * np.deg2rad(angle))
            positions = twist.apply(aligned) + target
            # Trial copies let us edit topology for scoring without altering
            # the original chain or the next trial's starting monomer.
            candidate = monomer.copy()
            for node, position in zip(nodes, positions):
                candidate.nodes[node].update(zip(('x', 'y', 'z'), map(float, position)))
            # Score the actual edited atoms: omit H_O4 and include O_a/H_Oa.
            candidate_chain = chain.copy()
            edit_beta_O_4_topology(candidate_chain, candidate, beta_atom, oxygen_atom)
            # Linker-specific exemptions: graph distances across the new bond
            # of at most three bonds (1-2, 1-3, and 1-4 pairs).
            chain_depth = nx.single_source_shortest_path_length(candidate_chain, beta_atom, cutoff=2)
            monomer_depth = nx.single_source_shortest_path_length(candidate, oxygen_atom, cutoff=2)
            # Each pair is (chain atom ID, incoming atom ID).
            # A depth is the number of existing bonds from the attachment atom.
            # Add 1 for the proposed C_b--O_4 bond between the two graphs.
            # Example: C_a -> C_b -> O_4 -> C_4 spans 1 + 1 + 1 = 3 bonds,
            # so that particular pair is exempt; other contacts still count.
            excluded_pairs = {
                (chain_atom, incoming_atom)
                for chain_atom, chain_distance in chain_depth.items()
                for incoming_atom, incoming_distance in monomer_depth.items()
                if chain_distance + 1 + incoming_distance <= 3
            }
            score = clash_score(candidate_chain, candidate, excluded_pairs=excluded_pairs)
            if score < best_score:
                best_score, best_positions = score, positions.copy()
            # Zero is the lowest possible score; no later rotation can improve it.
            if best_score == 0:
                break
        if best_score == 0:
            break

    # Commit only the winning coordinates. connect_beta_O_4 applies the actual
    # topology edits afterward; the edits above were only on scoring copies.
    for node, position in zip(nodes, best_positions):
        monomer.nodes[node].update(zip(('x', 'y', 'z'), map(float, position)))
    if best_score > 0:
        import warnings
        warnings.warn(
            "All beta-O-4 placements have close nonbonded contacts; "
            "using the least-overlapping candidate before optimization.",
            RuntimeWarning, stacklevel=2,
        )


def edit_beta_O_4_topology(chain, monomer, beta_atom, oxygen_atom):
    """Apply the specified beta-side and incoming O_4-side edits in place.

    Reduce C_a - C_b to order 1 and add O_a - H_Oa while retaining H_a.
    Seed the hydroxyl opposite the new C_b -> O_4 attachment direction.
    Attachment atoms are selected by connect_beta_O_4; template bonds are trusted.
    """
    # Restrict alpha-side edits to the same monomer as the selected beta carbon.
    monomer_id = chain.nodes[beta_atom]['monomer_id']
    alpha_atom = _find_attachment_atom(chain, 'C_a', monomer_id)
    
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
    # Keep H_a in place and point OH away from the newly attached monomer.
    hydroxyl_direction = (
        _coordinates(chain, beta_atom) - _coordinates(monomer, oxygen_atom)
    )
    direction_length = np.linalg.norm(hydroxyl_direction)
    if direction_length < 1e-10:
        raise ValueError("Cannot place hydroxyl: C_b and O_4 coordinates coincide.")
    hydroxyl_direction /= direction_length
    bond_distances = dictionaries.dictBondDistances()
    carbon_oxygen_distance = bond_distances[('C-O', 'alcohol')]
    oxygen_hydrogen_distance = bond_distances[('O-H', 'alcohol')]
    oxygen_position = alpha_position + carbon_oxygen_distance * hydroxyl_direction
    hydrogen_position = oxygen_position + oxygen_hydrogen_distance * hydroxyl_direction
    
    # Adjust bond order between alpha and beta carbons of "old" monomer
    chain.edges[alpha_atom, beta_atom]['bond_order'] = 1
    
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

