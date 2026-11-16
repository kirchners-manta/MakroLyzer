"""Check beta-1 chemistry and graph ownership independently of placement."""

import networkx as nx
import numpy as np
import pytest
from scipy.spatial.transform import Rotation

from MakroLyzer import dictionaries
from MakroLyzer.builder_modules import liginConnectionLogic as logic
from MakroLyzer.builder_modules.ligninBuilder import LigninBuilder


def monomer_pair(first='G', second='H'):
    builder = LigninBuilder.__new__(LigninBuilder)
    chain = builder.load_monomer(first, monomer_id=4, start_index=10)
    incoming = builder.load_monomer(second, monomer_id=9, start_index=200)
    # Start from distinct coordinates; placement must handle the translation.
    for node in incoming:
        incoming.nodes[node]['x'] += 3.0
    return builder, chain, incoming


@pytest.mark.parametrize('first', 'HGS')
@pytest.mark.parametrize('second', 'HGS')
def test_beta1_forms_closed_shell_spiro_ring_without_editing_inputs(first, second):
    builder, chain, incoming = monomer_pair(first, second)
    original_chain, original_incoming = chain.copy(), incoming.copy()
    original_ids = set(chain) | set(incoming)
    result = builder._connect_monomers(chain, 'beta-1', incoming, chain_monomer_id=4)

    assert nx.utils.graphs_equal(chain, original_chain)
    assert nx.utils.graphs_equal(incoming, original_incoming)
    assert nx.is_connected(result)
    assert len(result) == len(chain) + len(incoming) + 2
    for node, data in result.nodes(data=True):
        assert data['index'] == node
        assert data['degree'] == result.degree(node)
        valence = sum(result.edges[node, other]['bond_order'] for other in result[node])
        assert valence == {'C': 4, 'O': 2, 'H': 1}[data['element']]

    find = logic._find_attachment_atom
    bridge = find(result, 'O_a', 4)
    cycle = [find(result, 'C_1', 4), find(result, 'C_a', 4), bridge,
             find(result, 'C_a', 9), find(result, 'C_b', 9)]
    assert nx.is_isomorphic(result.subgraph(cycle), nx.cycle_graph(5))
    assert result.edges[find(result, 'C_4', 4), find(result, 'O_4', 4)]['bond_order'] == 2
    assert find(chain, 'H_O4', 4) not in result
    assert find(incoming, 'H_O4', 9) in result
    assert all(data['bond_order'] == 1 for _, _, data in result.edges(data=True)
               if data.get('linkage') == 'beta-1')
    for node in set(result) - original_ids:
        assert node > max(original_ids)
        assert result.nodes[node]['monomer_id'] == 4
        assert result.nodes[node]['monomer_type'] == first


@pytest.mark.parametrize('rotation_vector', [[0, 0, 0], [np.pi, 0, 0], [0.7, -1.2, 2]])
def test_beta1_seeds_beta_hydroxyl_along_a_normal_and_ether_at_alpha_midpoint(rotation_vector):
    _, chain, incoming = monomer_pair()
    rotation = Rotation.from_rotvec(rotation_vector)
    for graph in (chain, incoming):
        for node in graph:
            position = rotation.apply(logic._coordinates(graph, node)) + [5, -3, 7]
            graph.nodes[node].update(zip(('x', 'y', 'z'), position))
    find = logic._find_attachment_atom
    normal = logic.ring_normal(chain, 4)
    beta_position = logic._coordinates(chain, find(chain, 'C_b', 4))
    alpha_midpoint = (logic._coordinates(chain, find(chain, 'C_a', 4))
                      + logic._coordinates(incoming, find(incoming, 'C_a', 9))) / 2
    bridge = logic.edit_beta_1_topology(chain, incoming, find(chain, 'C_1', 4),
                                       find(incoming, 'C_b', 9))
    distances = dictionaries.dictBondDistances()
    oxygen = logic._coordinates(chain, find(chain, 'O_b', 4))
    hydrogen = logic._coordinates(chain, find(chain, 'H_Ob', 4))
    direction = logic.unit_vector(oxygen - beta_position)
    assert abs(np.dot(direction, normal)) == pytest.approx(1)
    assert np.linalg.norm(oxygen - beta_position) == pytest.approx(distances['C-O', 'alcohol'])
    np.testing.assert_allclose(hydrogen - oxygen, distances['O-H', 'alcohol'] * direction)
    np.testing.assert_allclose(logic._coordinates(chain, bridge), alpha_midpoint)
    # The helper must not create bare incoming nodes when forming the ether.
    assert not set(chain).intersection(incoming)
    assert set(chain[bridge]) == {find(chain, 'C_a', 4)}
    assert all(chain.nodes[node]['degree'] == degree for node, degree in chain.degree())


@pytest.mark.parametrize('preferred', [0, 1])
def test_beta1_selects_hydroxyl_side_using_complete_product_connectivity(monkeypatch, preferred):
    _, chain, incoming = monomer_pair('S', 'H')
    find = logic._find_attachment_atom
    c1, ca, cb = (find(chain, label, 4) for label in ('C_1', 'C_a', 'C_b'))
    incoming_ca, incoming_cb, incoming_ha = (
        find(incoming, label, 9) for label in ('C_a', 'C_b', 'H_a')
    )
    phenolic_hydrogen = find(chain, 'H_O4', 4)
    logic.place_beta_1(chain, incoming, c1, incoming_cb)
    incoming_positions = {node: logic._coordinates(incoming, node) for node in incoming}
    normal = logic.ring_normal(chain, 4)
    beta_position = logic._coordinates(chain, cb)
    trials = []

    def score(fixed, moving, excluded_pairs):
        oxygen = find(moving, 'O_b', 4)
        hydrogen = find(moving, 'H_Ob', 4)
        assert set(moving) == {oxygen, hydrogen}
        assert fixed.edges[c1, incoming_cb]['bond_order'] == 1
        assert fixed.edges[incoming_ca, incoming_cb]['bond_order'] == 1
        ether = find(fixed, 'O_a', 4)
        assert set(fixed[ether]) == {ca, incoming_ca}
        assert phenolic_hydrogen not in fixed
        # O_b--C_b--C_a--C_1 is exempt; the incoming alpha carbon
        # and its hydrogen are further away in the completed bond graph.
        assert (cb, oxygen) in excluded_pairs
        assert (c1, oxygen) in excluded_pairs
        assert (incoming_ca, oxygen) not in excluded_pairs
        assert (incoming_ha, oxygen) not in excluded_pairs
        trials.append(logic._coordinates(moving, oxygen))
        return 0.0 if len(trials) - 1 == preferred else 1.0

    monkeypatch.setattr(logic, 'clash_score', score)
    logic.edit_beta_1_topology(chain, incoming, c1, incoming_cb)
    assert len(trials) == 2
    np.testing.assert_allclose(trials[0] - beta_position, 1.43 * normal)
    np.testing.assert_allclose(trials[1] - beta_position, -1.43 * normal)
    np.testing.assert_allclose(logic._coordinates(chain, find(chain, 'O_b', 4)), trials[preferred])
    for node in incoming:
        np.testing.assert_array_equal(logic._coordinates(incoming, node), incoming_positions[node])


@pytest.mark.parametrize('invalid', ['chain_double_bond', 'incoming_double_bond', 'phenol', 'ring'])
def test_beta1_rejects_invalid_precursors_before_editing_either_graph(invalid):
    _, chain, incoming = monomer_pair()
    find = logic._find_attachment_atom
    if invalid.endswith('double_bond'):
        graph, identity = (chain, 4) if invalid.startswith('chain') else (incoming, 9)
        graph.edges[find(graph, 'C_a', identity), find(graph, 'C_b', identity)]['bond_order'] = 1
    elif invalid == 'phenol':
        chain.remove_node(find(chain, 'H_O4', 4))
    else:
        chain.remove_edge(find(chain, 'C_1', 4), find(chain, 'C_2', 4))
    original_chain, original_incoming = chain.copy(), incoming.copy()
    with pytest.raises(ValueError):
        logic.edit_beta_1_topology(chain, incoming, find(chain, 'C_1', 4),
                                   find(incoming, 'C_b', 9))
    assert nx.utils.graphs_equal(chain, original_chain)
    assert nx.utils.graphs_equal(incoming, original_incoming)


def test_beta1_only_edits_the_selected_chain_monomer():
    builder, chain, incoming = monomer_pair()
    spectator = builder.load_monomer('S', monomer_id=0, start_index=1000)
    chain = nx.compose(spectator, chain)
    result = logic.connect_beta_1(chain, incoming, 4)
    assert nx.utils.graphs_equal(result.subgraph(spectator), spectator)
    assert min(set(result) - set(chain) - set(incoming)) > max(chain)


@pytest.mark.parametrize('previous,following', [
    ('beta-1', '4-O-beta'), ('beta-1', 'beta-beta'), ('beta-1', 'beta-1'),
    ('beta-O-4', 'beta-1'), ('beta-5', 'beta-1'), ('beta-beta', 'beta-1'),
    ('4-O-beta', 'beta-1'),
])
def test_beta1_recipe_rejects_reusing_occupied_sites(previous, following):
    builder = LigninBuilder.__new__(LigninBuilder)
    builder.recipe = {'monomer': ['G', 'H', 'S'], 'linker': [previous, following]}
    with pytest.raises(ValueError, match='position of the shared monomer is already occupied'):
        builder.validate_recipe()


@pytest.mark.parametrize('following', ['beta-O-4', 'beta-5'])
def test_beta1_recipe_can_extend_using_free_ring_sites(following):
    builder = LigninBuilder.__new__(LigninBuilder)
    builder.recipe = {'monomer': ['G', 'H', 'S'], 'linker': ['beta-1', following]}
    builder.validate_recipe()
