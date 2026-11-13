"""Shared atom lookup, clash scoring, and geometry utilities for lignin builders."""

import networkx as nx
import numpy as np
from scipy.spatial import cKDTree
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
    # A missing or duplicated label is a template/selection error; choosing an
    # arbitrary match could silently connect the wrong atom. Return the node ID.
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
    # Without atoms on both sides there are no cross-graph contacts to score.
    if not fixed_graph or not moving_graph:
        return 0.0
    # Membership checks occur for every nearby pair; a set keeps them cheap.
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
    # Use a fresh float array in x/y/z order so callers can transform it
    # without changing the graph until coordinates are explicitly written back.
    position = np.array([graph.nodes[node][key] for key in ('x', 'y', 'z')], dtype=float)
    # Reject NaN and infinity before they reach rotations, SVD, or the KD-tree.
    if not np.isfinite(position).all():
        raise ValueError(f"Nonfinite coordinates for atom {node}.")
    return position


def ring_normal(graph, monomer_id):
    """Return the unit normal of the best-fit C_1 through C_6 ring plane.

    Select atoms by template label and monomer ID without modifying the graph.
    Orient the normal toward the cross product of C_1->C_2 and C_1->C_3.
    Raise ValueError for missing/ambiguous atoms, nonfinite coordinates, or a
    degenerate ring that cannot define a plane.
    """
    # Keep C_1 through C_6 in label order, independent of graph insertion order.
    ring = np.array([
        _coordinates(graph, _find_attachment_atom(graph, f'C_{i}', monomer_id))
        for i in range(1, 7)
    ])
    # Subtract the centroid so translation cannot affect the fitted plane.
    # SVD orders the principal spreads from largest to smallest. The final
    # row of axes is the unit direction with least spread: the plane normal.
    _, singular_values, axes = np.linalg.svd(ring - ring.mean(axis=0))
    # A plane needs two independent directions of spread. If the second
    # singular value vanishes, the points lie on a line or coincide.
    if singular_values[1] < 1e-10:
        raise ValueError("Cannot define a plane from a degenerate six-ring.")
    normal = axes[-1]
    # SVD leaves the sign arbitrary. Use the first three atoms to choose a
    # consistent side when their cross product has a nonzero normal component.
    reference = np.cross(ring[1] - ring[0], ring[2] - ring[0])
    if np.dot(normal, reference) < 0:
        normal = -normal
    return normal


def unit_vector(vector):
    """Return a normalized finite vector, rejecting degenerate directions."""
    # The Euclidean norm is the vector length; dividing by it keeps only
    # direction. Near-zero vectors have no reliable direction to normalize.
    vector = np.asarray(vector, dtype=float)
    length = np.linalg.norm(vector)
    if not np.isfinite(vector).all() or length < 1e-10:
        raise ValueError("Cannot normalize a nonfinite or degenerate direction.")
    return vector / length


def alignment_rotation(source, target):
    """Return a rotation mapping source onto target, including opposite vectors.

    Vector lengths do not matter. For opposite vectors choose a repeatable
    perpendicular axis; already aligned vectors need no rotation.
    """
    # For unit vectors, the cross-product length is sin(theta) and the dot
    # product is cos(theta). Their cross product gives the rotation axis.
    source, target = unit_vector(source), unit_vector(target)
    axis = np.cross(source, target)
    sine = np.linalg.norm(axis)
    # Roundoff can push a dot product slightly outside its mathematical range.
    cosine = np.clip(np.dot(source, target), -1.0, 1.0)
    # A rotation vector encodes unit axis times angle in radians. atan2
    # distinguishes acute and obtuse angles and remains stable near 180 degrees.
    if sine > 1e-10:
        return Rotation.from_rotvec(axis / sine * np.arctan2(sine, cosine))
    # Parallel and opposite vectors both have zero cross product. Opposite
    # vectors need a half-turn about any perpendicular axis. Choose the
    # Cartesian axis least aligned with source to avoid a tiny cross product.
    if cosine < 0:
        axis = unit_vector(np.cross(source, np.eye(3)[np.argmin(abs(source))]))
        return Rotation.from_rotvec(np.pi * axis)
    # Same direction: leave coordinates unchanged.
    return Rotation.identity()


def rotate_nodes(graph, nodes, rotation, pivot):
    """Rotate selected graph coordinates in place about a coordinate pivot."""
    # Rotation.apply rotates about the origin. Shift the pivot to the origin,
    # rotate, then shift back; the pivot itself therefore stays fixed.
    pivot = np.asarray(pivot, dtype=float)
    for node in nodes:
        position = rotation.apply(_coordinates(graph, node) - pivot) + pivot
        # Only coordinates change; atom attributes and graph bonds are retained.
        graph.nodes[node].update(zip(('x', 'y', 'z'), map(float, position)))


def branch_nodes(graph, left, right):
    """Return the right-side component after cutting a noncyclic bond.

    The graph is unchanged. Reject cyclic bonds because the right-side group
    cannot move independently of the left atom.
    """
    # Cut a temporary copy so selecting a movable group does not edit topology.
    cut = graph.copy()
    cut.remove_edge(left, right)
    # Connectivity collects the complete branch, including its hydrogens.
    nodes = nx.node_connected_component(cut, right)
    # If left is still reachable, an alternate path closes a cycle: there
    # is no independent right-side branch across this bond.
    if left in nodes:
        raise ValueError("Cannot isolate a branch across a cyclic bond.")
    return nodes


def ring_center(graph, monomer_id):
    """Return the centroid of the selected monomer's C_1 through C_6 atoms."""
    # Average the six carbon positions equally (axis=0 averages atoms, keeping
    # x/y/z). This is a geometric centroid, rather than a mass-weighted center.
    return np.mean([
        _coordinates(graph, _find_attachment_atom(graph, f'C_{i}', monomer_id))
        for i in range(1, 7)
    ], axis=0)
