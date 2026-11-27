"""Check directional 5-O-4 product construction through builder workflows."""

import json
from collections import Counter
from itertools import product

import networkx as nx
import numpy as np
import pytest

from MakroLyzer.builder_modules import liginConnectionLogic as logic
from MakroLyzer.builder_modules.ligninBuilder import LigninBuilder


def builder_for(tmp_path, monomers, linkers):
    """Write an ordinary JSON recipe and return its builder."""
    recipe = tmp_path / 'ether.json'
    recipe.write_text(json.dumps({'monomer': monomers, 'linker': linkers}))
    return LigninBuilder(recipe)


@pytest.mark.parametrize('chain_kind,incoming_kind', list(product('HGS', 'HG')))
def test_5_o_4_dimer_chemistry_and_relaxed_geometry(tmp_path, chain_kind, incoming_kind):
    graph = builder_for(tmp_path, [chain_kind, incoming_kind], ['biphenyl-ether-5-O-4']).run()
    atom = lambda label, identity: logic._find_attachment_atom(graph, label, identity)
    # Direct substitution retains the existing phenolic oxygen and both
    # monolignol side chains; coupling removes only phenolic H and C5-H.
    formulas = {'H': Counter(C=9, H=10, O=2), 'G': Counter(C=10, H=12, O=3),
                'S': Counter(C=11, H=14, O=4)}
    expected = formulas[chain_kind] + formulas[incoming_kind] - Counter(H=2)
    assert Counter(data['element'] for _, data in graph.nodes(data=True)) == expected
    assert nx.is_connected(graph)
    assert set(graph[atom('O_4', 0)]) == {atom('C_4', 0), atom('C_5', 1)}
    assert graph.edges[atom('O_4', 0), atom('C_5', 1)]['linkage'] == 'biphenyl-ether-5-O-4'
    assert atom('H_O4', 1) in graph
    assert not any(data['monomer_id'] == 0 and data['template_label'] == 'H_O4'
                   or data['monomer_id'] == 1 and data['template_label'] == 'H_5'
                   for _, data in graph.nodes(data=True))
    for identity in (0, 1):
        assert graph.edges[atom('C_a', identity), atom('C_b', identity)]['bond_order'] == 2
    for node, data in graph.nodes(data=True):
        assert data['index'] == node
        assert data['degree'] == graph.degree(node)
        assert sum(graph.edges[node, other]['bond_order'] for other in graph[node]) == {
            'C': 4, 'O': 2, 'H': 1,
        }[data['element']]
    # Check the whole product, not only the new bond or UFF convergence.
    for left, right, data in graph.edges(data=True):
        length = np.linalg.norm(logic._coordinates(graph, left) - logic._coordinates(graph, right))
        if 'H' in {graph.nodes[left]['element'], graph.nodes[right]['element']}:
            assert 0.85 < length < 1.2
        elif data['bond_order'] == 2:
            assert 1.1 < length < 1.5
        else:
            assert 1.25 < length < 1.75


@pytest.mark.parametrize('rotation_vector', [[0, 0, 0], [np.pi, 0, 0], [0.7, -1.2, 2]])
def test_placement_follows_attachment_axis_and_preserves_inputs(tmp_path, rotation_vector):
    builder = builder_for(tmp_path, ['S', 'G'], ['biphenyl-ether-5-O-4'])
    chain = builder.load_monomer('S')
    incoming = builder.load_monomer('G', monomer_id=1, start_index=max(chain) + 1)
    logic.rotate_nodes(incoming, list(incoming), logic.Rotation.from_rotvec(rotation_vector),
                       logic.ring_center(incoming, 1))
    originals = chain.copy(), incoming.copy()
    find = logic._find_attachment_atom
    c4, oxygen, c5 = find(chain, 'C_4', 0), find(chain, 'O_4', 0), find(incoming, 'C_5', 1)
    placed = incoming.copy()
    lengths = {(a, b): np.linalg.norm(logic._coordinates(incoming, a) - logic._coordinates(incoming, b))
               for a, b in incoming.edges}
    logic.place_biphenyl_ether_5_O_4(chain, placed, c5, c4, oxygen)
    direction = logic.unit_vector(logic._coordinates(chain, oxygen) - logic._coordinates(chain, c4))
    target = (logic._coordinates(chain, oxygen)
              + logic.dictionaries.dictBondDistances()[('C-O', 'ether')] * direction)
    np.testing.assert_allclose(logic._coordinates(placed, c5), target, atol=1e-10)
    for (a, b), length in lengths.items():
        assert np.linalg.norm(logic._coordinates(placed, a) - logic._coordinates(placed, b)) == pytest.approx(length)
    result = logic.connect_biphenyl_ether_5_O_4(chain, incoming, 0)
    assert nx.utils.graphs_equal(chain, originals[0])
    assert nx.utils.graphs_equal(incoming, originals[1])
    assert set(result) == (set(chain) | set(incoming)) - {
        find(chain, 'H_O4', 0), find(incoming, 'H_5', 1),
    }


@pytest.mark.parametrize('chain_kind', 'HGS')
def test_recipe_rejects_incoming_s_at_validation(tmp_path, chain_kind):
    with pytest.raises(ValueError, match='incoming.*S|S.*5 position'):
        builder_for(tmp_path, [chain_kind, 'S'], ['biphenyl-ether-5-O-4']).validate_recipe()


@pytest.mark.parametrize('linkers', [
    ['4-O-beta', 'biphenyl-ether-5-O-4'], ['biphenyl-ether-5-O-4', 'phenylcoumaran-beta-5'],
])
def test_recipe_rejects_occupied_middle_sites(tmp_path, linkers):
    with pytest.raises(ValueError, match='already occupied'):
        builder_for(tmp_path, ['G'] * 3, linkers).validate_recipe()


@pytest.mark.parametrize('linkers', [
    ['biphenyl-ether-5-O-4', 'biphenyl-ether-5-O-4'], ['beta-O-4', 'biphenyl-ether-5-O-4'], ['biphenyl-ether-5-O-4', 'beta-O-4'],
    ['biphenyl-ether-5-O-4', '4-O-beta'], ['biphenyl-ether-5-O-4', 'resinol-beta-beta'],
    ['phenylcoumaran-beta-5', 'biphenyl-ether-5-O-4'],
    ['biphenyl-ether-5-O-4', '1,2-diarylpropane-1,3-diol-beta-1'],
])
def test_5_o_4_preserves_existing_units_during_valid_growth(tmp_path, linkers):
    graph = builder_for(tmp_path, ['G'] * (len(linkers) + 1), linkers).run()
    assert nx.is_connected(graph)
    assert {data['monomer_id'] for _, data in graph.nodes(data=True)} == set(range(len(linkers) + 1))
    for node, data in graph.nodes(data=True):
        assert sum(graph.edges[node, other]['bond_order'] for other in graph[node]) == {
            'C': 4, 'O': 2, 'H': 1,
        }[data['element']]


def test_connector_rejects_occupied_oxygen_without_editing(tmp_path):
    builder = builder_for(tmp_path, ['G', 'G'], ['4-O-beta'])
    chain = builder.run()
    incoming = builder.load_monomer('G', monomer_id=2, start_index=max(chain) + 1)
    originals = chain.copy(), incoming.copy()
    with pytest.raises(ValueError, match='free chain phenolic'):
        logic.connect_biphenyl_ether_5_O_4(chain, incoming, 1)
    assert nx.utils.graphs_equal(chain, originals[0])
    assert nx.utils.graphs_equal(incoming, originals[1])


def test_5_o_4_cli_writes_xyz(tmp_path, monkeypatch):
    from MakroLyzer import MakroLyzer

    builder_for(tmp_path, ['G', 'G'], ['biphenyl-ether-5-O-4'])
    output = tmp_path / 'ether.xyz'
    monkeypatch.setattr('sys.argv', ['MakroLyzer', '-buildLignin', str(tmp_path / 'ether.json'),
                                   '--lignin-file', str(output)])
    MakroLyzer.main()
    rows = output.read_text().splitlines()
    # G/G diaryl ether: C20H22O6, from loss of two hydrogens only.
    assert int(rows[0]) == 48
    assert len(rows[2:]) == 48
