"""Lignin linkage construction and geometry hooks, independent of the builder."""

import numpy as np
from scipy.spatial.transform import Rotation
import networkx as nx
from MakroLyzer import dictionaries
from MakroLyzer.builder_modules.ligninHelpers import (
    _coordinates, _find_attachment_atom, alignment_rotation, branch_nodes,
    clash_score, ring_center, ring_normal, rotate_nodes, unit_vector,
)


# -----------------------------------------------------------------------------------------------
# 4-O-beta --------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------

def connect_4_O_beta(chain, monomer, chain_monomer_id):
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
    place_4_O_beta(chain, monomer, beta_atom, oxygen_atom)
    edit_4_O_beta_topology(chain, monomer, beta_atom, oxygen_atom)
    
    # Add nodes and edges from the incomming monomer to the chain
    result = chain.copy()
    result.add_nodes_from(monomer.nodes(data=True))
    result.add_edges_from(monomer.edges(data=True))
    # Add the new 4-O-beta bond connecting chain and new monomer
    result.add_edge(beta_atom, oxygen_atom, bond_order=1, linkage='4-O-beta')
    
    # update degree
    for node, degree in result.degree():
        result.nodes[node]['degree'] = degree
        
    return result

def place_4_O_beta(chain, monomer, beta_atom, oxygen_atom):
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
    
    normal = ring_normal(chain, chain_id)

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
        rotation = alignment_rotation(direction, outward)
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
            edit_4_O_beta_topology(candidate_chain, candidate, beta_atom, oxygen_atom)
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

    # Commit only the winning coordinates. connect_4_O_beta applies the actual
    # topology edits afterward; the edits above were only on scoring copies.
    for node, position in zip(nodes, best_positions):
        monomer.nodes[node].update(zip(('x', 'y', 'z'), map(float, position)))
    if best_score > 0:
        import warnings
        warnings.warn(
            "All 4-O-beta placements have close nonbonded contacts; "
            "using the least-overlapping candidate before optimization.",
            RuntimeWarning, stacklevel=2,
        )


def edit_beta_ether_topology(beta_graph, oxygen_graph, beta_atom, oxygen_atom):
    """Apply the beta-side and O_4-side edits in place.

    Reduce C_a - C_b to order 1 and add O_a - H_Oa while retaining H_a.
    Seed the hydroxyl opposite the new C_b -> O_4 attachment direction.
    The first graph contains C_b; the second contains O_4. Template bonds are trusted.
    """
    # Restrict alpha-side edits to the same monomer as the selected beta carbon.
    monomer_id = beta_graph.nodes[beta_atom]['monomer_id']
    alpha_atom = _find_attachment_atom(beta_graph, 'C_a', monomer_id)
    
    # Resolve the phenolic hydrogen within the selected O_4-side monomer.
    oxygen_monomer_id = oxygen_graph.nodes[oxygen_atom]['monomer_id']
    phenolic_hydrogen = _find_attachment_atom(oxygen_graph, 'H_O4', oxygen_monomer_id)
    
    # Require the expected starting bond before changing the topology.
    if (
        not beta_graph.has_edge(alpha_atom, beta_atom)
        or beta_graph.edges[alpha_atom, beta_atom].get('bond_order') != 2
    ):
        raise ValueError("Expected a C_a--C_b double bond before beta-ether edits.")

    # Get coordinates for OH, which will be added to alpha carbon of the beta-side monomer
    oxygen_id = max(max(beta_graph.nodes), max(oxygen_graph.nodes)) + 1
    hydrogen_id = oxygen_id + 1
    monomer_type = beta_graph.nodes[alpha_atom]['monomer_type']
    alpha_position = _coordinates(beta_graph, alpha_atom)
    # Keep H_a in place and point OH away from the O_4 side.
    hydroxyl_direction = (
        _coordinates(beta_graph, beta_atom) - _coordinates(oxygen_graph, oxygen_atom)
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
    
    # Adjust bond order between alpha and beta carbons of the beta-side monomer
    beta_graph.edges[alpha_atom, beta_atom]['bond_order'] = 1
    
    # Remove H_O4 from the O_4-side monomer
    oxygen_graph.remove_node(phenolic_hydrogen)
    
    # Add OH nodes and edges to the beta-side monomer
    for node, label, element, position in (
        (oxygen_id, 'O_a', 'O', oxygen_position),
        (hydrogen_id, 'H_Oa', 'H', hydrogen_position),
    ):
        beta_graph.add_node(
            node, index=node, element=element, atom_label=label,
            template_label=label, monomer_id=monomer_id,
            monomer_type=monomer_type,
            x=float(position[0]), y=float(position[1]), z=float(position[2]),
        )
    # Form the alpha hydroxyl group: C_a - O_a - H_Oa.
    beta_graph.add_edge(alpha_atom, oxygen_id, bond_order=1)
    beta_graph.add_edge(oxygen_id, hydrogen_id, bond_order=1)
    for node, degree in beta_graph.degree():
        beta_graph.nodes[node]['degree'] = degree
    for node, degree in oxygen_graph.degree():
        oxygen_graph.nodes[node]['degree'] = degree

def edit_4_O_beta_topology(chain, monomer, beta_atom, oxygen_atom):
    """Add alpha OH on the chain and remove H_O4 from the incoming monomer."""
    edit_beta_ether_topology(
        beta_graph=chain, oxygen_graph=monomer,
        beta_atom=beta_atom, oxygen_atom=oxygen_atom,
    )


# -----------------------------------------------------------------------------------------------
# beta-O-4 --------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------

def connect_beta_O_4(chain, monomer, chain_monomer_id):
    """Scaffold: chain O-4 connects to incoming monomer beta carbon.

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

    # beta belongs to the preceding imcomming monomer.
    # O_4 belongs to the recipe monomer. 
    # These variables hold node IDs, not atom data.
    oxygen_atom = _find_attachment_atom(chain, 'O_4', chain_monomer_id)
    beta_atom = _find_attachment_atom(monomer, 'C_b', incoming_id)
    
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
    Place the incoming monomer's C_b above the chain's O_4.

    Position C_b one ether bond length from O_4, opposite C_4,
    giving C_b--O_4--C_4 an initial angle of 180 degrees.
    Rigidly rotate the incoming monomer so its ring normal is
    parallel to the existing O_4--C_4 bond.

    Try both ring-normal signs and 20-degree rotations about the
    new bond. Choose the least-overlapping candidate, including
    the incoming monomer's new alpha hydroxyl.
    """
    chain_id = chain.nodes[oxygen_atom]['monomer_id']
    incoming_id = monomer.nodes[beta_atom]['monomer_id']
    
    normal = ring_normal(monomer, incoming_id)

    # Get the C4-O bond of the "old" monomer
    carbon = _find_attachment_atom(chain, 'C_4', chain_id)
    origin = _coordinates(chain, oxygen_atom)
    direction = _coordinates(chain, carbon) - origin
    length = np.linalg.norm(direction)
    
    direction /= length 

    bond_length = dictionaries.dictBondDistances()[('C-O', 'ether')]
    target = origin - bond_length * direction # target position of C_b

    nodes = list(monomer)
    beta_origin = _coordinates(monomer, beta_atom)
    original = np.array([_coordinates(monomer, node) for node in nodes]) - beta_origin

    best_score = float('inf')
    best_positions = None
    
    for sign in (1, -1):
        source = sign * normal
        rotation = alignment_rotation(source, direction)

        aligned = rotation.apply(original)

        for angle in range(0, 360, 20):
            twist = Rotation.from_rotvec(direction * np.deg2rad(angle))
            positions = twist.apply(aligned) + target

            candidate = monomer.copy()
            for node, position in zip(nodes, positions):
                candidate.nodes[node].update(
                    zip(('x', 'y', 'z'), map(float, position))
                )

            candidate_chain = chain.copy()
            edit_beta_O_4_topology(
                candidate_chain, candidate, beta_atom, oxygen_atom
            )

            # The chain attachment is O_4; the incoming attachment is C_b.
            chain_depth = nx.single_source_shortest_path_length(candidate_chain, oxygen_atom, cutoff=2)
            monomer_depth = nx.single_source_shortest_path_length(candidate, beta_atom, cutoff=2)
            excluded_pairs = {
                (chain_node, incoming_node)
                for chain_node, chain_distance in chain_depth.items()
                for incoming_node, incoming_distance in monomer_depth.items()
                if chain_distance + 1 + incoming_distance <= 3
            }

            score = clash_score(
                candidate_chain, candidate, excluded_pairs=excluded_pairs
            )
            if score < best_score:
                best_score = score
                best_positions = positions.copy()

            if best_score == 0:
                break
        if best_score == 0:
            break

    for node, position in zip(nodes, best_positions):
        monomer.nodes[node].update(
            zip(('x', 'y', 'z'), map(float, position))
        )

    if best_score > 0:
        import warnings
        warnings.warn(
            "All beta-O-4 placements have close nonbonded contacts; "
            "using the least-overlapping candidate before optimization.",
            RuntimeWarning, stacklevel=2,
        )
    
    
def edit_beta_O_4_topology(chain, monomer, beta_atom, oxygen_atom):
    """Add alpha OH on the incoming monomer and remove H_O4 from the chain."""
    edit_beta_ether_topology(
        beta_graph=monomer, oxygen_graph=chain,
        beta_atom=beta_atom, oxygen_atom=oxygen_atom,
    )

# -----------------------------------------------------------------------------------------------
# beta-5 ----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------

def connect_beta_5(chain, monomer, chain_monomer_id):
    """Scaffold: chain C5 connects to incoming monomer beta carbon.
                 chains C4 and incoming monomer C_a form an ether.
                 -> We keep the chain monomers O4 for the ether.

    Find corresponding atoms, place new monomer at the right position of the chain and 
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
    
    # imcomming monomer : C_a, C_5
    # chain monomer : C_5, O_C4 
    # These variables hold node IDs, not atom data.
    alpha_atom = _find_attachment_atom(monomer, 'C_a', incoming_id)
    beta_atom = _find_attachment_atom(monomer, 'C_b', incoming_id)
    C_5 = _find_attachment_atom(chain, 'C_5', chain_monomer_id)
    O_4 = _find_attachment_atom(chain, 'O_4', chain_monomer_id)
    
    # Place new monomer at the correct position near the existing chain
    place_beta_5(chain, monomer, C_5, O_4, alpha_atom, beta_atom)
    edit_beta_5_topology(chain, monomer, C_5, alpha_atom, beta_atom)
    
    # Add nodes and edges from the incomming monomer to the chain
    result = chain.copy()
    result.add_nodes_from(monomer.nodes(data=True))
    result.add_edges_from(monomer.edges(data=True))
    # Add the two new bonds connecting chain and new monomer
    result.add_edge(beta_atom, C_5, bond_order=1, linkage='beta-5')
    result.add_edge(alpha_atom, O_4, bond_order=1, linkage='beta-5')
    
    # update degree
    for node, degree in result.degree():
        result.nodes[node]['degree'] = degree
        
    return result
    
def place_beta_5(chain, monomer, C_5, O_4, alpha_atom, beta_atom):
    """
    Place incomming monomer in the same layer as the chain monomer.
    This is done by adjusting the ring normal of the incomming monomer such, that is is parallel
    to the one of the chain monomer.
    
    Make C_b-C_g bond of the incomming monomer parallel to its C_a-Cb bond to avoid clashes.
    Align the C_a-C_b bond of the incomming monomer parallel to C_4-C_5 bond of the chain monomer
    with a distance of 1.3*C-C bond lenth.
    Compare both ring-normal signs and retain the less-clashing orientation.
    """
    # get monomer IDs
    monomer_id = monomer.nodes[beta_atom]['monomer_id']
    chain_id = chain.nodes[C_5]['monomer_id']
    # get chain - C_4 and monomer C_g 
    carbon_4 = _find_attachment_atom(chain, 'C_4', chain_id)
    gamma_atom = _find_attachment_atom(monomer, 'C_g', monomer_id)

    # 1. Rotate the incoming ring normal onto the chain ring normal.
    normal = ring_normal(chain, chain_id)
    incoming_normal = ring_normal(monomer, monomer_id)
    rotate_nodes(
        monomer, list(monomer),
        alignment_rotation(incoming_normal, normal),
        ring_center(monomer, monomer_id),
    )

    # Keep the attachment bond in that plane. Templates can have a small
    # out-of-plane component; bend only the beta-side group to remove it.
    alpha = _coordinates(monomer, alpha_atom)
    direction = _coordinates(monomer, beta_atom) - alpha
    in_plane = direction - np.dot(direction, normal) * normal
    rotate_nodes(
        monomer, branch_nodes(monomer, alpha_atom, beta_atom),
        alignment_rotation(direction, in_plane), alpha,
    )

    # 2. Point C_b->C_g along C_a->C_b, keeping the gamma group rigid.
    # This deliberately changes the C_a--C_b--C_g angle to 180 degrees.
    beta = _coordinates(monomer, beta_atom)
    direction = beta - alpha
    rotate_nodes(
        monomer, branch_nodes(monomer, beta_atom, gamma_atom),
        alignment_rotation(_coordinates(monomer, gamma_atom) - beta, direction), beta,
    )

    # 3. Rotate within the common ring plane to align C_a->C_b with C_4->C_5.
    c4 = _coordinates(chain, carbon_4)
    c5 = _coordinates(chain, C_5)
    edge = c5 - c4
    edge_in_plane = unit_vector(edge - np.dot(edge, normal) * normal)
    direction = unit_vector(direction)
    angle = np.arctan2(np.dot(normal, np.cross(direction, edge_in_plane)),
                       np.dot(direction, edge_in_plane))
    rotate_nodes(monomer, list(monomer), Rotation.from_rotvec(angle * normal), alpha)

    # Translate outside the C_4--C_5 edge at 1.3 times its C--C length.
    # The ring normals are already aligned, so both planes share this normal.
    # Their centers locate each plane along that normal (its "height").
    incoming_center = ring_center(monomer, monomer_id)
    chain_center = ring_center(chain, chain_id)

    # Start halfway between the incoming attachment atoms C_a and C_b.
    midpoint = (_coordinates(monomer, alpha_atom) + _coordinates(monomer, beta_atom)) / 2
    # dot(midpoint - center, normal) is the signed height above the ring plane.
    # Subtract that height along the unit normal to project onto the plane.
    # This projected midpoint is the incoming reference point we will move.
    anchor = midpoint - np.dot(midpoint - incoming_center, normal) * normal

    # Likewise, project the chain's C_4--C_5 midpoint onto its ring plane.
    # Both reference points now lie in their respective ring planes, so moving
    # anchor onto a point in the chain plane makes the parallel planes coincide.
    target = (c4 + c5) / 2
    target -= np.dot(target - chain_center, normal) * normal

    # Crossing the plane normal with the edge direction gives a direction
    # within the plane, perpendicular to C_4--C_5. Normalize it to unit length.
    outward = unit_vector(np.cross(normal, edge_in_plane))
    # There are two perpendicular directions. Choose the one pointing away
    # from the chain ring center toward the edge, rather than into the ring.
    if np.dot(outward, target - chain_center) < 0:
        outward = -outward

    # The desired anchor position is beyond the edge by 1.3 * its bond length.
    # Destination minus current anchor gives the required translation vector.
    # This offset sets the separation of the projected bond midpoints; it does
    # not set either of the two new attachment bond lengths directly.
    shift = target + 1.3 * np.linalg.norm(edge) * outward - anchor
    # Apply the same shift to every incoming atom: internal geometry and
    # orientation stay unchanged, and the existing chain stays fixed.
    for node in monomer:
        position = _coordinates(monomer, node) + shift
        monomer.nodes[node].update(zip(('x', 'y', 'z'), map(float, position)))

    # Score the product atoms using copies; remove hydrogens only on these copies.
    scoring_chain = chain.copy()
    scoring_monomer = monomer.copy()
    edit_beta_5_topology(scoring_chain, scoring_monomer, C_5, alpha_atom, beta_atom)
    # Include both new bonds when calculating graph distances.
    joined = nx.compose(scoring_chain, scoring_monomer)
    joined.add_edge(beta_atom, C_5, bond_order=1)
    joined.add_edge(alpha_atom, O_4, bond_order=1)
    # Skip cross-monomer pairs within three bonds: their proximity is expected.
    excluded_pairs = {
        (fixed, moving)
        for fixed in scoring_chain
        for moving in nx.single_source_shortest_path_length(joined, fixed, cutoff=3)
        if moving in scoring_monomer
    }
    # Preserve the current orientation as the first candidate.
    nodes = list(monomer)
    original = np.array([_coordinates(monomer, node) for node in nodes])
    # Flip the incoming ring to the other side of the C_a--C_b axis.
    pivot = _coordinates(monomer, alpha_atom)
    axis = unit_vector(_coordinates(monomer, beta_atom) - pivot)
    flipped = Rotation.from_rotvec(np.pi * axis).apply(original - pivot) + pivot
    # Small template deviations can put the axis slightly above the ring
    # plane. Restore the common plane height after flipping.
    flipped_center = Rotation.from_rotvec(np.pi * axis).apply(
        ring_center(monomer, monomer_id) - pivot
    ) + pivot
    flipped -= np.dot(flipped_center - chain_center, normal) * normal
    # Compare both orientations and keep the one with less nonbonded overlap.
    best_score = float('inf')
    best_positions = original
    for positions in (original, flipped):
        candidate = scoring_monomer.copy()
        for node, position in zip(nodes, positions):
            candidate.nodes[node].update(zip(('x', 'y', 'z'), map(float, position)))
        score = clash_score(scoring_chain, candidate, excluded_pairs=excluded_pairs)
        if score < best_score:
            best_score, best_positions = score, positions
    # Commit only the winning coordinates; actual topology edits happen later.
    for node, position in zip(nodes, best_positions):
        monomer.nodes[node].update(zip(('x', 'y', 'z'), map(float, position)))


def edit_beta_5_topology(chain, monomer, C_5, alpha_atom, beta_atom):
    """
    incoming monomer: reduce C_a - C_b bond order
    chain: remove H5 (H or G), remove H_O4
    """
    # get monomer IDs
    monomer_id = monomer.nodes[beta_atom]['monomer_id']
    chain_monomer_id = chain.nodes[C_5]['monomer_id']
    
    ### incoming monomer ###
    # C_a-C_b bond order
    if (
        not monomer.has_edge(alpha_atom, beta_atom)
        or monomer.edges[alpha_atom, beta_atom].get('bond_order') != 2
    ):
        raise ValueError("Expected a C_a--C_b double bond before beta-ether edits.")
    monomer.edges[alpha_atom, beta_atom]['bond_order'] = 1
    
    ### chain monomer ###
    # remove H5
    H_5 = _find_attachment_atom(chain, 'H_5', chain_monomer_id)
    chain.remove_node(H_5)
    # H_O4
    H_O4 = _find_attachment_atom(chain, 'H_O4', chain_monomer_id)
    chain.remove_node(H_O4)
    
    for node, degree in chain.degree():
        chain.nodes[node]['degree'] = degree
    for node, degree in monomer.degree():
        monomer.nodes[node]['degree'] = degree