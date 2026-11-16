"""Lignin linkage construction and geometry hooks, independent of the builder."""

import numpy as np
from scipy.spatial.transform import Rotation
import networkx as nx
from MakroLyzer import dictionaries
from MakroLyzer.builder_modules import globalHelpers
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
        
# -----------------------------------------------------------------------------------------------
# beta-beta -------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------

def connect_beta_beta(chain, monomer, chain_monomer_id, *, max_iterations=1500):
    """Scaffold: chain C_b connects to incoming monomer C_b.
                 chains O_g connects to incoming monomer C_a.
                 Incoming monomers O_g connects to chains C_a.

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
    
    # imcomming monomer : C_b, C_a, O_g
    # chain monomer : C_b, C_a, O_g
    # These variables hold node IDs, not atom data.
    monomer_C_a = _find_attachment_atom(monomer, 'C_a', incoming_id)
    monomer_C_b = _find_attachment_atom(monomer, 'C_b', incoming_id)
    monomer_O_g = _find_attachment_atom(monomer, 'O_g', incoming_id)
    chain_C_a = _find_attachment_atom(chain, 'C_a', chain_monomer_id)
    chain_C_b = _find_attachment_atom(chain, 'C_b', chain_monomer_id)
    chain_O_g = _find_attachment_atom(chain, 'O_g', chain_monomer_id)
    
    # Place new monomer at the correct position near the existing chain
    place_beta_beta(chain, monomer, chain_C_a, chain_C_b, chain_O_g, monomer_C_a, monomer_C_b, monomer_O_g,
                    max_iterations=max_iterations)
    edit_beta_beta_topology(chain, monomer, chain_C_a, chain_C_b, monomer_C_a, monomer_C_b)
    
    # Add nodes and edges from the incomming monomer to the chain
    result = chain.copy()
    result.add_nodes_from(monomer.nodes(data=True))
    result.add_edges_from(monomer.edges(data=True))
    # Add the three new bonds connecting chain and new monomer
    result.add_edge(monomer_C_b, chain_C_b, bond_order=1, linkage='beta-beta')
    result.add_edge(monomer_C_a, chain_O_g, bond_order=1, linkage='beta-beta')
    result.add_edge(chain_C_a, monomer_O_g, bond_order=1, linkage='beta-beta')
    
    # update degree
    for node, degree in result.degree():
        result.nodes[node]['degree'] = degree
        
    return result

def place_beta_beta(chain, monomer, chain_C_a, chain_C_b, chain_O_g, monomer_C_a, monomer_C_b, monomer_O_g,
                    *, max_iterations=1500):
    """
    Place incoming monomer in the same layer as the chain monomer, similar to beta-5.
    Make C_1-C_3 vactors of incoming and chain monomer parallel.
    Make C_a-C_b-C_g linear and parallel to C_1-C_3, in the ring layer. -> For incoming monomer and 
    chain monomer individually. 
    Make C_b's lay on one line perpendicular to C1-C3 vectors with a distance of C-C "single" to each other.
    Rank eight orientations by relaxed product energy (two UFF passes per trial).
    Commit the starting coordinates of the winning trial; relaxation copies
    do not change the requested initial alignment.
    """
    # Select the two monomers and use the existing ring normal as reference.
    chain_id = chain.nodes[chain_C_a]['monomer_id']
    monomer_id = monomer.nodes[monomer_C_a]['monomer_id']
    normal = ring_normal(chain, chain_id)

    def ring_direction(graph, identity, ring_axis):
        c1 = _coordinates(graph, _find_attachment_atom(graph, 'C_1', identity))
        c3 = _coordinates(graph, _find_attachment_atom(graph, 'C_3', identity))
        # Project C_1->C_3 onto the fitted ring plane: real templates can
        # deviate slightly from a perfectly flat ring. Keep only its direction.
        direction = c3 - c1
        return unit_vector(direction - np.dot(direction, ring_axis) * ring_axis)

    def straighten(graph, identity, alpha, beta, direction):
        # Bend the complete downstream groups, retaining their bond lengths.
        # Moving only C_b or C_g would leave their hydrogens behind.
        gamma = _find_attachment_atom(graph, 'C_g', identity)
        origin = _coordinates(graph, alpha)
        # Parallel includes opposite directions: point away from this ring.
        if np.dot(direction, origin - ring_center(graph, identity)) < 0:
            direction = -direction
        # Pivot at C_a and rotate the complete beta-side group onto the axis.
        rotate_nodes(graph, branch_nodes(graph, alpha, beta),
                     alignment_rotation(_coordinates(graph, beta) - origin, direction), origin)
        # Pivot at the new C_b position and align the gamma-side group too.
        # C_a->C_b and C_b->C_g now point along the same line.
        origin = _coordinates(graph, beta)
        rotate_nodes(graph, branch_nodes(graph, beta, gamma),
                     alignment_rotation(_coordinates(graph, gamma) - origin, direction), origin)

    # 3: make C_a--C_b--C_g linear and parallel to C_1--C_3
    # for both monomers. Do this locally before aligning the two rings.
    # Straighten each side chain along its own ring's projected C_1->C_3 axis.
    along = ring_direction(chain, chain_id, normal)
    straighten(chain, chain_id, chain_C_a, chain_C_b, along)
    incoming_normal = ring_normal(monomer, monomer_id)
    incoming_along = ring_direction(monomer, monomer_id, incoming_normal)
    straighten(monomer, monomer_id, monomer_C_a, monomer_C_b, incoming_along)

    # Preparation of steps 1 and 2: each frame contains an in-plane axis,
    # its perpendicular in-plane axis, and the ring normal. Matching frames
    # therefore aligns both the ring planes and the C_1--C_3 directions.
    source_frame = np.column_stack((incoming_along,
                                    np.cross(incoming_normal, incoming_along), incoming_normal))
    nodes = list(monomer)
    # Put incoming C_b at the origin so trial rotations leave it fixed.
    beta_origin = _coordinates(monomer, monomer_C_b)
    original = np.array([_coordinates(monomer, node) for node in nodes]) - beta_origin
    chain_beta = _coordinates(chain, chain_C_b)
    # 4 direction: within the ring plane, perpendicular to C_1--C_3.
    across = unit_vector(np.cross(normal, along))
    distances = dictionaries.dictBondDistances()
    cc_length = distances[('C-C', 'single')]

    # Use the product topology for exemptions through any of the three bonds.
    # Remove both H_Og atoms and reduce the double bonds on temporary copies.
    scoring_chain, scoring_monomer = chain.copy(), monomer.copy()
    edit_beta_beta_topology(scoring_chain, scoring_monomer,
                            chain_C_a, chain_C_b, monomer_C_a, monomer_C_b)
    joined = nx.compose(scoring_chain, scoring_monomer)
    for left, right in ((chain_C_b, monomer_C_b), (chain_O_g, monomer_C_a),
                        (chain_C_a, monomer_O_g)):
        joined.add_edge(left, right, bond_order=1)
    # Pairs separated by up to three bonds are intentionally close; skip
    # them in the nonbonded clash score, including paths through ether bonds.
    excluded_pairs = {
        (fixed, moving)
        for fixed in scoring_chain
        for moving in nx.single_source_shortest_path_length(joined, fixed, cutoff=3)
        if moving in scoring_monomer
    }
    best_score, best_positions = float('inf'), None
    trials = []
    # Parallel axes may point in either direction. Test both, both normal
    # signs, and both sides of the beta atom to avoid overlapping the rings.
    for direction_sign in (1, -1):
        for normal_sign in (1, -1):
            target_along, target_normal = direction_sign * along, normal_sign * normal
            target_frame = np.column_stack((target_along,
                                            np.cross(target_normal, target_along), target_normal))
            # Protocol steps 1 and 2: rotate the entire incoming monomer.
            # Normals and C_1--C_3 axes become parallel (either sign is allowed).
            rotated = original @ (target_frame @ source_frame.T).T
            for side in (1, -1):
                # Protocol step 4: put incoming C_b one single C--C length
                # from chain C_b along the perpendicular in-plane direction.
                # This brings the ring layers together; small template heights
                # above their fitted planes can leave a slight plane offset.
                positions = rotated + chain_beta + side * cc_length * across
                candidate = scoring_monomer.copy()
                for node, position in zip(nodes, positions):
                    if node not in candidate:  # H_Og is absent from the product.
                        continue
                    candidate.nodes[node].update(zip(('x', 'y', 'z'), map(float, position)))
                # Clash scoring cannot assess strained bonded contacts, which
                # are excluded above. Keep its ranking, then compare the actual
                # relaxed products rather than tuning arbitrary closure weights.
                score = clash_score(scoring_chain, candidate, excluded_pairs=excluded_pairs)
                trials.append((score, positions.copy(), candidate))

    # Additional selection after the protocol: compare all eight complete
    # product geometries, because a clash score omits bonded strain.
    # Evaluate both normal signs even when one has a slightly worse clash score.
    # A short ether closure can otherwise trap UFF in a distorted local minimum.
    for _, positions, candidate in sorted(trials, key=lambda trial: trial[0]):
        product = nx.compose(scoring_chain, candidate)
        for left, right in ((chain_C_b, monomer_C_b), (chain_O_g, monomer_C_a),
                            (chain_C_a, monomer_O_g)):
            product.add_edge(left, right, bond_order=1)
        relaxed = globalHelpers.optimize_structure(product, max_iterations=max_iterations)
        # A second relaxation can escape a stalled first pass. Compare the
        # lower energy of the two passes for this candidate.
        retried = globalHelpers.optimize_structure(relaxed, max_iterations=max_iterations)
        score = min(relaxed.graph['uff_optimization']['final_energy'],
                    retried.graph['uff_optimization']['final_energy'])
        if score < best_score:
            # Commit the initial placement, preserving the requested straight
            # side chains. The builder performs its normal optimization later.
            best_score, best_positions = score, positions
    # Write back the winning initial placement, not the relaxed trial geometry.
    for node, position in zip(nodes, best_positions):
        monomer.nodes[node].update(zip(('x', 'y', 'z'), map(float, position)))


def edit_beta_beta_topology(chain, monomer, chain_C_a, chain_C_b, monomer_C_a, monomer_C_b):
    """
    incoming monomer: reduce C_a - C_b bond order, remove H_Og
    chain monomer: : reduce C_a - C_b bond order, remove H_Og
    """
    # get monomer IDs
    monomer_id = monomer.nodes[monomer_C_a]['monomer_id']
    chain_monomer_id = chain.nodes[chain_C_a]['monomer_id']
    
    ### incoming monomer ###
    # C_a-C_b bond order
    if (
        not monomer.has_edge(monomer_C_a, monomer_C_b)
        or monomer.edges[monomer_C_a, monomer_C_b].get('bond_order') != 2
    ):
        raise ValueError("Expected a C_a--C_b double bond before beta-ether edits.")
    monomer.edges[monomer_C_a, monomer_C_b]['bond_order'] = 1
    # remove H_Og
    H_Og = _find_attachment_atom(monomer, 'H_Og', monomer_id)
    monomer.remove_node(H_Og)
    
    ### chain monomer ###
    # C_a-C_b bond order
    if (
        not chain.has_edge(chain_C_a, chain_C_b)
        or chain.edges[chain_C_a, chain_C_b].get('bond_order') != 2
    ):
        raise ValueError("Expected a C_a--C_b double bond before beta-ether edits.")
    chain.edges[chain_C_a, chain_C_b]['bond_order'] = 1
    # remove H_Og
    H_Og = _find_attachment_atom(chain, 'H_Og', chain_monomer_id)
    chain.remove_node(H_Og)
    
    for node, degree in chain.degree():
        chain.nodes[node]['degree'] = degree
    for node, degree in monomer.degree():
        monomer.nodes[node]['degree'] = degree


# -----------------------------------------------------------------------------------------------
# beta-1 ----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------

def connect_beta_1(chain, monomer, chain_monomer_id):
    """
    Connect incoming C_b to chain C_1 and close the C_a--O--C_a ether.

    Return a new graph without editing the inputs.
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
    if any(data['monomer_id'] == incoming_id for _, data in chain.nodes(data=True)):
        raise ValueError("Incoming monomer ID must not already occur in the chain.")

    # Incoming monomer: C_a, C_b; chain monomer: C_1.
    # These variables hold node IDs, not atom data.
    monomer_C_a = _find_attachment_atom(monomer, 'C_a', incoming_id)
    monomer_C_b = _find_attachment_atom(monomer, 'C_b', incoming_id)
    chain_C_1 = _find_attachment_atom(chain, 'C_1', chain_monomer_id)

    # Place new monomer at the correct position near the existing chain
    place_beta_1(chain, monomer, chain_C_1, monomer_C_b)
    ether_oxygen = edit_beta_1_topology(chain, monomer, chain_C_1, monomer_C_b)

    # Add nodes and edges from the incoming monomer to the chain
    result = chain.copy()
    result.add_nodes_from(monomer.nodes(data=True))
    result.add_edges_from(monomer.edges(data=True))
    # Add the two new bonds connecting chain and new monomer
    result.add_edge(monomer_C_b, chain_C_1, bond_order=1, linkage='beta-1')
    result.add_edge(monomer_C_a, ether_oxygen, bond_order=1, linkage='beta-1')

    # update degree
    for node, degree in result.degree():
        result.nodes[node]['degree'] = degree

    return result

def place_beta_1(chain, monomer, chain_C_1, monomer_C_b):
    """
    Place the incoming monomer above the chain with parallel ring planes.
    Align the incoming ring normal with the chain ring normal.

    For the incoming and the chain monomer, we need to rearrange the carbon tails at C1.
    For this we mirror the positions of the carbon atoms of the rings at the mirror plane spanned
    by the C1-C6 vector and the ring normal vector.
    Ca comes on C2', Cb on C3', Cg on C4'. Preserve the Cg-Og bond length
    and point Cg->Og parallel to C2->C1 in each monomer.
    Then we arrange the incoming monomer such that its Ca-Cb vector is parallel to the chain's Ca-C1 vector,
    align the projected midpoints of these vectors, and separate the midpoints
    by 1.3*C-C "single" along the positive chain ring normal. Both atom pairs
    are vertically aligned for matching planar attachment edges. Preserve tiny
    endpoint mismatches from unequal edge lengths or nonplanar XYZ templates.
    """
    # get monomer IDs
    monomer_id = monomer.nodes[monomer_C_b]['monomer_id']
    chain_id = chain.nodes[chain_C_1]['monomer_id']

    # 1. Rotate the incoming ring normal onto the chain ring normal.
    normal = ring_normal(chain, chain_id)
    incoming_normal = ring_normal(monomer, monomer_id)
    rotate_nodes(
        monomer, list(monomer),
        alignment_rotation(incoming_normal, normal),
        ring_center(monomer, monomer_id),
    )

    # 2. Reflect C_2..C_4 across the plane through C_1 spanned by
    # C_1->C_6 and the ring normal, and put C_a, C_b, C_g there.
    def mirror_tail(graph, identity):
        c1 = _find_attachment_atom(graph, 'C_1', identity)
        c6 = _find_attachment_atom(graph, 'C_6', identity)
        origin = _coordinates(graph, c1)
        mirror_normal = unit_vector(np.cross(
            _coordinates(graph, c6) - origin, ring_normal(graph, identity),
        ))
        targets = []
        for i in range(2, 5):
            position = _coordinates(graph, _find_attachment_atom(graph, f'C_{i}', identity))
            targets.append(position - 2 * np.dot(position - origin, mirror_normal) * mirror_normal)

        backbone = [c1] + [_find_attachment_atom(graph, label, identity)
                           for label in ('C_a', 'C_b', 'C_g', 'O_g')]
        downstream = [branch_nodes(graph, left, right)
                      for left, right in zip(backbone, backbone[1:])]
        original = {node: _coordinates(graph, node) for node in downstream[0] | {c1}}
        gamma, oxygen = backbone[-2:]
        co_length = np.linalg.norm(original[oxygen] - original[gamma])
        c2 = _find_attachment_atom(graph, 'C_2', identity)
        targets.append(targets[-1] + co_length * unit_vector(origin - _coordinates(graph, c2)))
        positions = [origin] + targets
        for i, right in enumerate(backbone[1:], 1):
            # Move each atom with its own substituents. Fit the two local
            # backbone directions so hydrogens follow the bend instead of
            # pointing into the next carbon. Backbone anchors stay exact.
            left = backbone[i - 1]
            source = [unit_vector(original[left] - original[right])]
            target = [unit_vector(positions[i - 1] - positions[i])]
            nodes = downstream[i - 1]
            if i < len(backbone) - 1:
                source.append(unit_vector(original[backbone[i + 1]] - original[right]))
                target.append(unit_vector(positions[i + 1] - positions[i]))
                rotation, _ = Rotation.align_vectors(target, source)
                nodes = nodes - downstream[i]
            else:
                rotation = alignment_rotation(source[0], target[0])
            for node in nodes:
                position = rotation.apply(original[node] - original[right]) + positions[i]
                graph.nodes[node].update(zip(('x', 'y', 'z'), map(float, position)))

    mirror_tail(chain, chain_id)
    mirror_tail(monomer, monomer_id)

    # 3. Rotate about the common ring normal to align incoming C_a->C_b
    # with chain C_a->C_1. Use projections because XYZ rings have small
    # departures from planarity; keep the exact mirrored tail positions.
    chain_C_a = _find_attachment_atom(chain, 'C_a', chain_id)
    monomer_C_a = _find_attachment_atom(monomer, 'C_a', monomer_id)
    chain_alpha = _coordinates(chain, chain_C_a)
    incoming_alpha = _coordinates(monomer, monomer_C_a)
    direction = _coordinates(monomer, monomer_C_b) - incoming_alpha
    edge = _coordinates(chain, chain_C_1) - chain_alpha
    along = unit_vector(edge - np.dot(edge, normal) * normal)
    incoming_along = unit_vector(direction - np.dot(direction, normal) * normal)
    angle = np.arctan2(np.dot(normal, np.cross(incoming_along, along)),
                       np.dot(incoming_along, along))
    rotate_nodes(monomer, list(monomer), Rotation.from_rotvec(angle * normal), incoming_alpha)

    # 4. Stack the incoming monomer above the chain along the ring normal.
    # Center the attachment edges in projection rather than displacing them
    # sideways. Unequal template edge lengths leave small endpoint residuals;
    # the attachment midpoints are separated by exactly distance * normal.
    distance = 1.3 * dictionaries.dictBondDistances()[('C-C', 'single')]
    chain_midpoint = (chain_alpha + _coordinates(chain, chain_C_1)) / 2
    incoming_midpoint = (incoming_alpha + _coordinates(monomer, monomer_C_b)) / 2
    shift = chain_midpoint + distance * normal - incoming_midpoint
    for node in monomer:
        position = _coordinates(monomer, node) + shift
        monomer.nodes[node].update(zip(('x', 'y', 'z'), map(float, position)))

def edit_beta_1_topology(chain, monomer, chain_C_1, monomer_C_b):
    """Apply the beta-1 scaffold edits in place and return the new ether O ID.

    Reduce both C_a--C_b double bonds. On the chain side, form the dienone
    ring (C_2=C_3, C_5=C_6 and C_4=O_4), remove H_O4, and add O_b--H_Ob
    at C_b along the less-clashing ring-normal direction using the alcohol
    bond distances. Score both hydroxyl orientations against the completed
    product, excluding atom pairs separated by up to three bonds.

    Store the new O_a at the midpoint of the two alpha carbons in ``chain``
    and bond it to chain C_a. ``connect_beta_1`` adds its bond to incoming
    C_a and the C_1--C_b bond after merging the two graphs. This avoids
    creating incomplete copies of incoming atoms in the chain graph.
    """
    monomer_id = monomer.nodes[monomer_C_b]['monomer_id']
    chain_monomer_id = chain.nodes[chain_C_1]['monomer_id']
    monomer_C_a = _find_attachment_atom(monomer, 'C_a', monomer_id)
    chain_C_a = _find_attachment_atom(chain, 'C_a', chain_monomer_id)
    chain_C_b = _find_attachment_atom(chain, 'C_b', chain_monomer_id)
    # Resolve and validate both sides before making any edits.
    for graph, alpha, beta in (
        (monomer, monomer_C_a, monomer_C_b),
        (chain, chain_C_a, chain_C_b),
    ):
        if not graph.has_edge(alpha, beta) or graph.edges[alpha, beta].get('bond_order') != 2:
            raise ValueError("Expected a C_a--C_b double bond before beta-1 edits.")

    H_O4 = _find_attachment_atom(chain, 'H_O4', chain_monomer_id)
    chain_C_2 = _find_attachment_atom(chain, 'C_2', chain_monomer_id)
    chain_C_3 = _find_attachment_atom(chain, 'C_3', chain_monomer_id)
    chain_C_4 = _find_attachment_atom(chain, 'C_4', chain_monomer_id)
    chain_C_5 = _find_attachment_atom(chain, 'C_5', chain_monomer_id)
    chain_C_6 = _find_attachment_atom(chain, 'C_6', chain_monomer_id)
    chain_O_4 = _find_attachment_atom(chain, 'O_4', chain_monomer_id)
    ring_bonds = (
        (chain_C_1, chain_C_2, 1), (chain_C_2, chain_C_3, 2),
        (chain_C_3, chain_C_4, 1), (chain_C_4, chain_C_5, 1),
        (chain_C_5, chain_C_6, 2), (chain_C_6, chain_C_1, 1),
        (chain_C_4, chain_O_4, 2),
    )
    if any(not chain.has_edge(left, right) for left, right, _ in ring_bonds):
        raise ValueError("Expected an intact six-ring and C_4--O_4 bond before beta-1 edits.")

    # Allocate across both graphs before deleting H_O4, so IDs cannot collide
    # with incoming atoms or reuse the ID of the removed hydrogen.
    oxygen_id = max(max(chain), max(monomer)) + 1
    hydrogen_id = oxygen_id + 1
    ether_oxygen_id = oxygen_id + 2
    monomer_type = chain.nodes[chain_C_a]['monomer_type']
    distances = dictionaries.dictBondDistances()
    normal = ring_normal(chain, chain_monomer_id)
    oxygen_position = _coordinates(chain, chain_C_b) + distances[('C-O', 'alcohol')] * normal
    hydrogen_position = oxygen_position + distances[('O-H', 'alcohol')] * normal
    ether_position = (_coordinates(chain, chain_C_a) + _coordinates(monomer, monomer_C_a)) / 2

    monomer.edges[monomer_C_a, monomer_C_b]['bond_order'] = 1
    chain.edges[chain_C_a, chain_C_b]['bond_order'] = 1
    chain.remove_node(H_O4)
    for left, right, order in ring_bonds:
        chain.edges[left, right]['bond_order'] = order

    # Both the beta hydroxyl and the alpha ether belong to the selected chain monomer.
    for node, label, element, position in (
        (oxygen_id, 'O_b', 'O', oxygen_position),
        (hydrogen_id, 'H_Ob', 'H', hydrogen_position),
        (ether_oxygen_id, 'O_a', 'O', ether_position),
    ):
        chain.add_node(
            node, index=node, element=element, atom_label=label,
            template_label=label, monomer_id=chain_monomer_id,
            monomer_type=monomer_type,
            x=float(position[0]), y=float(position[1]), z=float(position[2]),
        )
    chain.add_edge(chain_C_b, oxygen_id, bond_order=1)
    chain.add_edge(oxygen_id, hydrogen_id, bond_order=1)
    chain.add_edge(chain_C_a, ether_oxygen_id, bond_order=1, linkage='beta-1')

    # Include the actual cross-monomer bonds when exempting nearby bonded
    # atoms. Score only the movable hydroxyl against the rest of the product,
    # so other close contacts cannot obscure the difference between sides.
    product = nx.compose(chain, monomer)
    product.add_edge(chain_C_1, monomer_C_b, bond_order=1)
    product.add_edge(ether_oxygen_id, monomer_C_a, bond_order=1)
    hydroxyl_nodes = {oxygen_id, hydrogen_id}
    fixed_nodes = set(product) - hydroxyl_nodes
    excluded_pairs = {
        (fixed, moving)
        for moving in hydroxyl_nodes
        for fixed in nx.single_source_shortest_path_length(product, moving, cutoff=3)
        if fixed in fixed_nodes
    }
    beta_position = _coordinates(chain, chain_C_b)
    offsets = (oxygen_position - beta_position, hydrogen_position - beta_position)
    best_score, best_positions = float('inf'), None
    for sign in (1, -1):
        positions = [beta_position + sign * offset for offset in offsets]
        for node, position in zip((oxygen_id, hydrogen_id), positions):
            product.nodes[node].update(zip(('x', 'y', 'z'), map(float, position)))
        score = clash_score(product.subgraph(fixed_nodes), product.subgraph(hydroxyl_nodes),
                            excluded_pairs=excluded_pairs)
        if score < best_score:
            best_score, best_positions = score, positions
    for node, position in zip((oxygen_id, hydrogen_id), best_positions):
        chain.nodes[node].update(zip(('x', 'y', 'z'), map(float, position)))

    for node, degree in chain.degree():
        chain.nodes[node]['degree'] = degree
    for node, degree in monomer.degree():
        monomer.nodes[node]['degree'] = degree

    return ether_oxygen_id
