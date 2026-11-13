"""Verify the beta-5 alignment, side-chain bending, and offset protocol."""

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from MakroLyzer.builder_modules import liginConnectionLogic as logic
from MakroLyzer.builder_modules.ligninBuilder import LigninBuilder


@pytest.mark.parametrize('rotation_vector', [[0, 0, 0], [np.pi, 0, 0], [0.7, -1.2, 2]])
def test_beta5_placement_follows_protocol(rotation_vector):
    builder = LigninBuilder.__new__(LigninBuilder)
    chain = builder.load_monomer('G')
    monomer = builder.load_monomer('G', monomer_id=1, start_index=max(chain) + 1)
    rotation = Rotation.from_rotvec(rotation_vector)
    for node in monomer:
        position = rotation.apply(logic._coordinates(monomer, node)) + [3, -5, 7]
        monomer.nodes[node].update(zip(('x', 'y', 'z'), position))

    def positions(graph):
        return np.array([logic._coordinates(graph, n) for n in graph])

    def atom(graph, label, monomer_id):
        return logic._find_attachment_atom(graph, label, monomer_id)

    before = positions(monomer)
    chain_before = positions(chain)
    alpha, beta = (atom(monomer, label, 1) for label in ('C_a', 'C_b'))
    c4, c5, o4 = (atom(chain, label, 0) for label in ('C_4', 'C_5', 'O_4'))
    logic.place_beta_5(chain, monomer, c5, o4, alpha, beta)
    after = positions(monomer)
    np.testing.assert_array_equal(positions(chain), chain_before)
    # Bending changes internal angles, but every existing bond length stays fixed.
    indices = {node: i for i, node in enumerate(monomer)}
    for left, right in monomer.edges:
        i, j = indices[left], indices[right]
        assert np.linalg.norm(after[i] - after[j]) == pytest.approx(
            np.linalg.norm(before[i] - before[j]), abs=1e-12
        )
    normal = logic.ring_normal(chain, 0)
    assert abs(np.dot(logic.ring_normal(monomer, 1), normal)) == pytest.approx(1)
    ring_in = np.array([logic._coordinates(monomer, atom(monomer, f'C_{i}', 1)) for i in range(1, 7)])
    ring_chain = np.array([logic._coordinates(chain, atom(chain, f'C_{i}', 0)) for i in range(1, 7)])
    assert np.dot(ring_in.mean(axis=0) - ring_chain.mean(axis=0), normal) == pytest.approx(0, abs=1e-12)
    incoming = logic._coordinates(monomer, beta) - logic._coordinates(monomer, alpha)
    edge = logic._coordinates(chain, c5) - logic._coordinates(chain, c4)
    projected = [v - np.dot(v, normal) * normal for v in (incoming, edge)]
    np.testing.assert_allclose(projected[0] / np.linalg.norm(projected[0]), projected[1] / np.linalg.norm(projected[1]), atol=1e-12)
    gamma = atom(monomer, 'C_g', 1)
    gamma_direction = logic._coordinates(monomer, gamma) - logic._coordinates(monomer, beta)
    np.testing.assert_allclose(
        gamma_direction / np.linalg.norm(gamma_direction),
        incoming / np.linalg.norm(incoming), atol=1e-12,
    )
    assert np.dot(incoming, normal) == pytest.approx(0, abs=1e-12)
    midpoint = (logic._coordinates(monomer, alpha) + logic._coordinates(monomer, beta)) / 2
    chain_midpoint = (logic._coordinates(chain, c4) + logic._coordinates(chain, c5)) / 2
    separation = midpoint - chain_midpoint
    separation -= np.dot(separation, normal) * normal
    assert np.linalg.norm(separation) == pytest.approx(1.3 * np.linalg.norm(edge))
    assert np.dot(separation, projected[1]) == pytest.approx(0, abs=1e-12)
    assert np.dot(separation, chain_midpoint - ring_chain.mean(axis=0)) > 0


@pytest.mark.parametrize('preferred', [0, 1])
def test_beta5_selects_lower_score_using_edited_topology(monkeypatch, preferred):
    builder = LigninBuilder.__new__(LigninBuilder)
    chain = builder.load_monomer('G')
    monomer = builder.load_monomer('H', monomer_id=1, start_index=max(chain) + 1)
    find = logic._find_attachment_atom
    alpha, beta = (find(monomer, label, 1) for label in ('C_a', 'C_b'))
    c5, o4 = (find(chain, label, 0) for label in ('C_5', 'O_4'))
    trials = []

    def score(fixed, moving, excluded_pairs):
        assert find(chain, 'H_5', 0) not in fixed
        assert find(chain, 'H_O4', 0) not in fixed
        assert moving.edges[alpha, beta]['bond_order'] == 1
        assert {(c5, beta), (o4, alpha), (c5, alpha), (o4, beta)} <= excluded_pairs
        trials.append(np.array([logic._coordinates(moving, n) for n in monomer]))
        return 0.0 if len(trials) - 1 == preferred else 1.0

    monkeypatch.setattr(logic, 'clash_score', score)
    logic.place_beta_5(chain, monomer, c5, o4, alpha, beta)
    assert len(trials) == 2
    np.testing.assert_allclose(
        [logic._coordinates(monomer, n) for n in monomer], trials[preferred], atol=1e-12,
    )
    assert not np.allclose(trials[0], trials[1])
    # Trial topology edits must not leak into either original graph.
    assert find(chain, 'H_5', 0) in chain
    assert find(chain, 'H_O4', 0) in chain
    assert monomer.edges[alpha, beta]['bond_order'] == 2


def test_beta5_dimer_optimizes_to_reasonable_attachment_lengths():
    from MakroLyzer.builder_modules import globalHelpers

    builder = LigninBuilder.__new__(LigninBuilder)
    chain = builder.load_monomer('G')
    monomer = builder.load_monomer('H', monomer_id=1, start_index=max(chain) + 1)
    result = globalHelpers.optimize_structure(logic.connect_beta_5(chain, monomer, 0))
    for left, left_id, right, right_id, expected in (
        ('C_b', 1, 'C_5', 0, 1.5), ('C_a', 1, 'O_4', 0, 1.42),
    ):
        a = logic._find_attachment_atom(result, left, left_id)
        b = logic._find_attachment_atom(result, right, right_id)
        distance = np.linalg.norm(logic._coordinates(result, a) - logic._coordinates(result, b))
        assert distance == pytest.approx(expected, abs=0.1)
    diagnostics = result.graph['uff_optimization']
    assert diagnostics['final_energy'] < diagnostics['initial_energy']
