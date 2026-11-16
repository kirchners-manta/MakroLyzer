"""Verify the reflected-tail beta-1 placement protocol before relaxation."""

import json

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from MakroLyzer import dictionaries
from MakroLyzer.builder_modules import liginConnectionLogic as logic
from MakroLyzer.builder_modules.ligninBuilder import LigninBuilder


@pytest.mark.parametrize('first', 'HGS')
@pytest.mark.parametrize('second', 'HGS')
@pytest.mark.parametrize('rotation_vector', [[0, 0, 0], [np.pi, 0, 0], [0.7, -1.2, 2]])
def test_beta1_placement_follows_reflection_alignment_and_distance_protocol(first, second, rotation_vector):
    builder = LigninBuilder.__new__(LigninBuilder)
    chain = builder.load_monomer(first, monomer_id=4, start_index=10)
    incoming = builder.load_monomer(second, monomer_id=9, start_index=200)
    spectator = builder.load_monomer('S', monomer_id=0, start_index=500)
    chain.add_nodes_from(spectator.nodes(data=True))
    chain.add_edges_from(spectator.edges(data=True))
    for graph, rotation, translation in (
        (chain, Rotation.from_rotvec([0.4, -0.8, 0.2]), [5, 2, -3]),
        (incoming, Rotation.from_rotvec(rotation_vector), [-2, 7, 11]),
    ):
        for node in graph:
            position = rotation.apply(logic._coordinates(graph, node)) + translation
            graph.nodes[node].update(zip(('x', 'y', 'z'), position))
    originals = [chain.copy(), incoming.copy()]
    find = logic._find_attachment_atom
    c1, ca = (find(chain, label, 4) for label in ('C_1', 'C_a'))
    incoming_ca, incoming_cb = (find(incoming, label, 9) for label in ('C_a', 'C_b'))
    tail = logic.branch_nodes(chain, c1, ca)
    logic.place_beta_1(chain, incoming, c1, incoming_cb)

    # The existing ring, its substituents, and other monomers stay fixed.
    for node in set(chain) - tail:
        np.testing.assert_array_equal(logic._coordinates(chain, node), logic._coordinates(originals[0], node))
    for graph, original, identity in zip((chain, incoming), originals, (4, 9)):
        assert dict(graph.edges) == dict(original.edges)
        for node in graph:
            assert {k: v for k, v in graph.nodes[node].items() if k not in ('x', 'y', 'z')} == {
                k: v for k, v in original.nodes[node].items() if k not in ('x', 'y', 'z')}
        p = lambda label: logic._coordinates(graph, find(graph, label, identity))
        origin = p('C_1')
        mirror_normal = logic.unit_vector(np.cross(p('C_6') - origin, logic.ring_normal(graph, identity)))
        reflection = np.eye(3) - 2 * np.outer(mirror_normal, mirror_normal)
        for label, index in zip(('C_a', 'C_b', 'C_g'), range(2, 5)):
            expected = origin + reflection @ (p(f'C_{index}') - origin)
            np.testing.assert_allclose(p(label), expected, atol=1e-12)
        np.testing.assert_allclose(logic.unit_vector(p('O_g') - p('C_g')),
                                   logic.unit_vector(p('C_1') - p('C_2')), atol=1e-12)
        # Mapping onto ring positions changes three backbone lengths.
        # C_g--O_g, all C--H and O--H bonds retain their original lengths.
        backbone = [find(graph, label, identity) for label in ('C_1', 'C_a', 'C_b', 'C_g')]
        moved_bonds = {frozenset(pair) for pair in zip(backbone, backbone[1:])}
        for left, right in graph.edges:
            if frozenset((left, right)) not in moved_bonds:
                before = logic._coordinates(original, left) - logic._coordinates(original, right)
                after = logic._coordinates(graph, left) - logic._coordinates(graph, right)
                assert np.linalg.norm(after) == pytest.approx(np.linalg.norm(before), abs=1e-12)

    normal = logic.ring_normal(chain, 4)
    np.testing.assert_allclose(logic.ring_normal(incoming, 9), normal, atol=1e-12)
    expected_distance = 1.3 * dictionaries.dictBondDistances()['C-C', 'single']
    assert np.dot(logic.ring_center(incoming, 9) - logic.ring_center(chain, 4), normal) == pytest.approx(expected_distance, abs=1e-4)
    directions = [logic._coordinates(chain, c1) - logic._coordinates(chain, ca),
                  logic._coordinates(incoming, incoming_cb) - logic._coordinates(incoming, incoming_ca)]
    projected = [logic.unit_vector(v - np.dot(v, normal) * normal) for v in directions]
    np.testing.assert_allclose(projected[0], projected[1], atol=1e-12)
    # Tiny template nonplanarity remains in the exact mirrored coordinates.
    assert np.dot(*(logic.unit_vector(v) for v in directions)) > 1 - 1e-8
    displacements = []
    for fixed, moving in ((c1, incoming_cb), (ca, incoming_ca)):
        displacement = logic._coordinates(incoming, moving) - logic._coordinates(chain, fixed)
        displacements.append(displacement)
        assert np.dot(displacement, normal) == pytest.approx(expected_distance, abs=1e-4)
        lateral = displacement - np.dot(displacement, normal) * normal
        assert np.linalg.norm(lateral) < 0.01
    np.testing.assert_allclose(np.mean(displacements, axis=0), expected_distance * normal, atol=1e-12)


def test_beta1_handles_equal_length_edges_in_planar_regular_rings():
    builder = LigninBuilder.__new__(LigninBuilder)
    chain = builder.load_monomer('H')
    incoming = builder.load_monomer('G', monomer_id=1, start_index=100)
    find = logic._find_attachment_atom
    for graph, identity in ((chain, 0), (incoming, 1)):
        for i in range(1, 7):
            theta = (i - 1) * np.pi / 3
            graph.nodes[find(graph, f'C_{i}', identity)].update(
                x=1.4 * np.cos(theta), y=1.4 * np.sin(theta), z=0.0,
            )
    logic.place_beta_1(chain, incoming, find(chain, 'C_1', 0), find(incoming, 'C_b', 1))
    np.testing.assert_allclose(logic.ring_normal(chain, 0), logic.ring_normal(incoming, 1), atol=1e-12)
    for left, right in (('C_1', 'C_b'), ('C_a', 'C_a')):
        vector = (logic._coordinates(incoming, find(incoming, right, 1))
                  - logic._coordinates(chain, find(chain, left, 0)))
        np.testing.assert_allclose(vector, 1.3 * 1.54 * logic.ring_normal(chain, 0), atol=1e-12)


@pytest.mark.parametrize('first', 'HGS')
@pytest.mark.parametrize('second', 'HGS')
def test_beta1_dimer_relaxes_to_reasonable_attachment_lengths(first, second):
    from MakroLyzer.builder_modules import globalHelpers

    builder = LigninBuilder.__new__(LigninBuilder)
    chain = builder.load_monomer(first)
    incoming = builder.load_monomer(second, monomer_id=1, start_index=100)
    result = globalHelpers.optimize_structure(logic.connect_beta_1(chain, incoming, 0))
    for left, right, data in result.edges(data=True):
        if data.get('linkage') == 'beta-1':
            expected = 1.54 if result.nodes[left]['element'] == result.nodes[right]['element'] == 'C' else 1.42
            distance = np.linalg.norm(logic._coordinates(result, left) - logic._coordinates(result, right))
            assert distance == pytest.approx(expected, abs=0.15)
    assert result.graph['uff_optimization']['final_energy'] < result.graph['uff_optimization']['initial_energy']


def test_sh_recipe_returns_relaxed_bonds_through_full_builder(tmp_path):
    recipe = tmp_path / 'sh.json'
    recipe.write_text(json.dumps({'monomer': ['S', 'H'], 'linker': ['beta-1']}))
    result = LigninBuilder(recipe).run()
    find = logic._find_attachment_atom
    for left, left_id, right, right_id, expected in (
        ('C_1', 0, 'C_b', 1, 1.54),
        ('C_a', 0, 'O_a', 0, 1.42),
        ('O_a', 0, 'C_a', 1, 1.42),
    ):
        a, b = find(result, left, left_id), find(result, right, right_id)
        length = np.linalg.norm(logic._coordinates(result, a) - logic._coordinates(result, b))
        assert length == pytest.approx(expected, abs=0.1)
    # Check the whole returned structure, beyond just the three new bonds.
    for a, b, data in result.edges(data=True):
        elements = {result.nodes[a]['element'], result.nodes[b]['element']}
        length = np.linalg.norm(logic._coordinates(result, a) - logic._coordinates(result, b))
        if 'H' in elements:
            assert 0.85 < length < 1.2
        elif data['bond_order'] == 2:
            assert 1.1 < length < 1.5
        else:
            assert 1.25 < length < 1.75
