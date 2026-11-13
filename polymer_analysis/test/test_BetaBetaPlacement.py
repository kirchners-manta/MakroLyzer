"""Check beta-beta placement geometry and the three optimized closures."""

import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from MakroLyzer.builder_modules import globalHelpers, liginConnectionLogic as logic
from MakroLyzer.builder_modules.ligninBuilder import LigninBuilder


@pytest.mark.parametrize('rotation_vector', [[0, 0, 0], [np.pi, 0, 0], [0.7, -1.2, 2]])
def test_beta_beta_placement_protocol(rotation_vector):
    builder = LigninBuilder.__new__(LigninBuilder)
    chain = builder.load_monomer('G')
    incoming = builder.load_monomer('H', monomer_id=1, start_index=max(chain) + 1)
    rotation = Rotation.from_rotvec(rotation_vector)
    for node in incoming:
        position = rotation.apply(logic._coordinates(incoming, node)) + [2, -3, 5]
        incoming.nodes[node].update(zip(('x', 'y', 'z'), position))
    original_chain = {n: dict(d) for n, d in chain.nodes(data=True)}
    lengths = {(id(g), u, v): np.linalg.norm(logic._coordinates(g, u) - logic._coordinates(g, v))
               for g in (chain, incoming) for u, v in g.edges}
    result = logic.connect_beta_beta(chain, incoming, 0)
    assert dict(chain.nodes(data=True)) == original_chain

    def p(label, identity):
        return logic._coordinates(result, logic._find_attachment_atom(result, label, identity))

    normal = logic.ring_normal(result, 0)
    assert abs(np.dot(normal, logic.ring_normal(result, 1))) == pytest.approx(1)
    ring_axis = p('C_3', 0) - p('C_1', 0)
    ring_axis -= np.dot(ring_axis, normal) * normal
    ring_axis = logic.unit_vector(ring_axis)
    delta = p('C_b', 1) - p('C_b', 0)
    assert np.linalg.norm(delta) == pytest.approx(1.54)
    assert np.dot(delta, normal) == pytest.approx(0, abs=1e-12)
    assert np.dot(delta, ring_axis) == pytest.approx(0, abs=1e-12)
    for identity in (0, 1):
        ab = logic.unit_vector(p('C_b', identity) - p('C_a', identity))
        bg = logic.unit_vector(p('C_g', identity) - p('C_b', identity))
        np.testing.assert_allclose(ab, bg, atol=1e-12)
        assert abs(np.dot(ab, ring_axis)) == pytest.approx(1)
    for g in (chain, incoming):
        for u, v in g.edges:
            if u in result and v in result:
                assert np.linalg.norm(logic._coordinates(result, u) - logic._coordinates(result, v)) == pytest.approx(lengths[id(g), u, v])


@pytest.mark.parametrize('first', 'HGS')
@pytest.mark.parametrize('second', 'HGS')
def test_beta_beta_optimized_closures(first, second):
    builder = LigninBuilder.__new__(LigninBuilder)
    chain = builder.load_monomer(first)
    incoming = builder.load_monomer(second, monomer_id=1, start_index=max(chain) + 1)
    result = globalHelpers.optimize_structure(logic.connect_beta_beta(chain, incoming, 0))
    relaxed = globalHelpers.optimize_structure(result)
    if relaxed.graph['uff_optimization']['final_energy'] < result.graph['uff_optimization']['final_energy']:
        result = relaxed
    for u, v, data in result.edges(data=True):
        if data.get('linkage') == 'beta-beta':
            expected = 1.54 if result.nodes[u]['element'] == result.nodes[v]['element'] == 'C' else 1.42
            assert np.linalg.norm(logic._coordinates(result, u) - logic._coordinates(result, v)) == pytest.approx(expected, abs=0.15)
