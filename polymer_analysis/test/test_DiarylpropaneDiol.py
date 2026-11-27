"""Exercise direct acyclic beta-1 construction through reachable recipes."""

import json
from collections import Counter
from itertools import product

import networkx as nx
import numpy as np
import pytest

from MakroLyzer.builder_modules import liginConnectionLogic as logic
from MakroLyzer.builder_modules.ligninBuilder import LigninBuilder


DIOL = '1,2-diarylpropane-1,3-diol-beta-1'


def builder_for(tmp_path, monomers, linkers):
    """Return a builder using the same JSON contract as the CLI."""
    recipe = tmp_path / 'diol.json'
    recipe.write_text(json.dumps({'monomer': monomers, 'linker': linkers}))
    return LigninBuilder(recipe)


@pytest.mark.parametrize('kinds', list(product('HGS', repeat=2)))
def test_diol_formula_connectivity_and_relaxed_geometry(tmp_path, kinds):
    graph = builder_for(tmp_path, list(kinds), [DIOL]).run()
    atom = lambda label, identity: logic._find_attachment_atom(graph, label, identity)
    # The retained product loses the first monolignol's C3H5O tail and
    # gains OH on incoming alpha. Ring substitution preserves aromaticity.
    formulas = {'H': Counter(C=9, H=10, O=2), 'G': Counter(C=10, H=12, O=3),
                'S': Counter(C=11, H=14, O=4)}
    expected = formulas[kinds[0]] + formulas[kinds[1]] - Counter(C=3, H=4)
    assert Counter(data['element'] for _, data in graph.nodes(data=True)) == expected
    assert nx.is_connected(graph)
    assert set(graph[atom('C_1', 0)]) == {atom('C_2', 0), atom('C_6', 0), atom('C_b', 1)}
    assert set(graph[atom('O_a', 1)]) == {atom('C_a', 1), atom('H_Oa', 1)}
    assert graph.edges[atom('C_a', 1), atom('C_b', 1)]['bond_order'] == 1
    for node, data in graph.nodes(data=True):
        assert data['index'] == node
        assert data['degree'] == graph.degree(node)
        assert sum(graph.edges[node, other]['bond_order'] for other in graph[node]) == {
            'C': 4, 'O': 2, 'H': 1,
        }[data['element']]
    for left, right, data in graph.edges(data=True):
        length = np.linalg.norm(logic._coordinates(graph, left) - logic._coordinates(graph, right))
        if 'H' in {graph.nodes[left]['element'], graph.nodes[right]['element']}:
            assert 0.85 < length < 1.2
        elif data['bond_order'] == 2:
            assert 1.1 < length < 1.5
        else:
            assert 1.25 < length < 1.75


def test_placement_follows_requested_vectors_without_changing_chain(tmp_path):
    builder = builder_for(tmp_path, ['G', 'S'], [DIOL])
    chain = builder.load_monomer('G')
    incoming = builder.load_monomer('S', monomer_id=1, start_index=max(chain) + 1)
    original = chain.copy()
    find = logic._find_attachment_atom
    c1 = find(chain, 'C_1', 0)
    beta = find(incoming, 'C_b', 1)
    gamma = find(incoming, 'C_g', 1)
    # Give the incoming template an arbitrary orientation, as after a
    # previous placement; the normal alignment must not assume XYZ axes.
    logic.rotate_nodes(incoming, list(incoming), logic.Rotation.from_euler('xyz', [40, 70, 15], degrees=True),
                       logic.ring_center(incoming, 1))
    gamma_branch = logic.branch_nodes(incoming, beta, gamma)
    lengths = {node: np.linalg.norm(logic._coordinates(incoming, node)
                                   - logic._coordinates(incoming, gamma)) for node in gamma_branch}
    logic.place_diarylpropane_diol_beta_1(chain, incoming, c1, beta)
    assert nx.utils.graphs_equal(chain, original)
    np.testing.assert_allclose(logic.ring_normal(chain, 0), logic.ring_normal(incoming, 1), atol=1e-10)
    np.testing.assert_allclose(logic._coordinates(incoming, beta),
                               logic._coordinates(chain, find(chain, 'C_a', 0)), atol=1e-10)
    np.testing.assert_allclose(
        logic.unit_vector(logic._coordinates(incoming, gamma) - logic._coordinates(incoming, beta)),
        logic.unit_vector(logic._coordinates(incoming, find(incoming, 'C_a', 1))
                          - logic._coordinates(incoming, find(incoming, 'C_1', 1))), atol=1e-10,
    )
    for node, length in lengths.items():
        assert np.linalg.norm(logic._coordinates(incoming, node)
                              - logic._coordinates(incoming, gamma)) == pytest.approx(length)


@pytest.mark.parametrize('linkers', [[DIOL, 'beta-O-4'], [DIOL, 'phenylcoumaran-beta-5'],
                                   ['4-O-beta', DIOL], ['4-O-beta', DIOL, 'beta-O-4']])
def test_diol_preserves_previous_units_and_allows_ring_growth(tmp_path, linkers):
    graph = builder_for(tmp_path, ['G'] * (len(linkers) + 1), linkers).run()
    assert nx.is_connected(graph)
    assert {data['monomer_id'] for _, data in graph.nodes(data=True)} == set(range(len(linkers) + 1))


@pytest.mark.parametrize('preceding', ['beta-O-4', 'phenylcoumaran-beta-5', 'resinol-beta-beta', DIOL])
def test_recipe_rejects_removing_attached_units(tmp_path, preceding):
    builder = builder_for(tmp_path, ['G'] * 3, [preceding, DIOL])
    with pytest.raises(ValueError, match='tail removal'):
        builder.validate_recipe()


@pytest.mark.parametrize('following', ['4-O-beta', 'resinol-beta-beta'])
def test_recipe_rejects_reusing_diol_beta(tmp_path, following):
    with pytest.raises(ValueError, match='beta position'):
        builder_for(tmp_path, ['G'] * 3, [DIOL, following]).validate_recipe()


def test_connector_rejects_tail_loss_without_modifying_inputs(tmp_path):
    builder = builder_for(tmp_path, ['G', 'H'], ['beta-O-4'])
    chain = builder.run()
    incoming = builder.load_monomer('S', monomer_id=2, start_index=max(chain) + 1)
    originals = chain.copy(), incoming.copy()
    with pytest.raises(ValueError, match='other lignin units'):
        logic.connect_diarylpropane_diol_beta_1(chain, incoming, 1)
    assert nx.utils.graphs_equal(chain, originals[0])
    assert nx.utils.graphs_equal(incoming, originals[1])


def test_connector_preserves_inputs_and_retained_ring_attributes(tmp_path):
    builder = builder_for(tmp_path, ['G', 'H'], ['4-O-beta'])
    chain = builder.run()
    incoming = builder.load_monomer('S', monomer_id=2, start_index=max(chain) + 1)
    originals = chain.copy(), incoming.copy()
    c1 = logic._find_attachment_atom(chain, 'C_1', 1)
    alpha = logic._find_attachment_atom(chain, 'C_a', 1)
    tail = logic.branch_nodes(chain, c1, alpha)
    graph = logic.connect_diarylpropane_diol_beta_1(chain, incoming, 1)
    assert nx.utils.graphs_equal(chain, originals[0])
    assert nx.utils.graphs_equal(incoming, originals[1])
    assert not set(graph).intersection(tail)
    for node in set(chain) - tail:
        # Connectivity-derived degree is refreshed; identities and coordinates
        # stay intact until the caller performs its normal UFF optimization.
        assert {key: value for key, value in graph.nodes[node].items() if key != 'degree'} == {
            key: value for key, value in chain.nodes[node].items() if key != 'degree'
        }


def test_diol_cli_writes_retained_product(tmp_path, monkeypatch):
    from MakroLyzer import MakroLyzer

    builder_for(tmp_path, ['G', 'G'], [DIOL])
    output = tmp_path / 'diol.xyz'
    monkeypatch.setattr('sys.argv', ['MakroLyzer', '-buildLignin', str(tmp_path / 'diol.json'),
                                   '--lignin-file', str(output)])
    MakroLyzer.main()
    rows = output.read_text().splitlines()
    # G/G retained diol is C17H20O6, rather than two complete C9 units.
    assert int(rows[0]) == 43
    assert len(rows[2:]) == 43
