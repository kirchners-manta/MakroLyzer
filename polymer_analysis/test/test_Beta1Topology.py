"""Check the native spirodienone trimer through real builder workflows."""

import json
from collections import Counter
from itertools import product

import networkx as nx
import pytest

from MakroLyzer.builder_modules import liginConnectionLogic as logic
from MakroLyzer.builder_modules.ligninBuilder import LigninBuilder


SPIRO = 'spirodienone-beta-1'


def build(tmp_path, monomers, linkers):
    """Write a user recipe and return its builder and completed structure."""
    recipe = tmp_path / 'recipe.json'
    recipe.write_text(json.dumps({'monomer': monomers, 'linker': linkers}))
    builder = LigninBuilder(recipe)
    return builder, builder.run()


@pytest.mark.parametrize('monomers', list(product('HGS', repeat=3)))
def test_spirodienone_trimer_connectivity_and_formula(tmp_path, monomers):
    _, graph = build(tmp_path, monomers, ['beta-O-4', SPIRO])
    find = logic._find_attachment_atom
    atom = lambda label, identity: find(graph, label, identity)

    # Elder et al., Figure 1: incoming unit 2 forms beta--1 and alpha--O--alpha
    # with unit 1, which retains its beta--O_4 connection to unit 0. These
    # tests use complete H/G/S units, including unit 0's C3 side chain.
    cycle = [atom('C_1', 1), atom('C_a', 1), atom('O_a', 1),
             atom('C_a', 2), atom('C_b', 2)]
    assert nx.is_isomorphic(graph.subgraph(cycle), nx.cycle_graph(5))
    assert set(graph[atom('O_a', 1)]) == {atom('C_a', 1), atom('C_a', 2)}
    assert set(graph[atom('C_1', 1)]) == {
        atom('C_2', 1), atom('C_6', 1), atom('C_a', 1), atom('C_b', 2),
    }
    assert set(graph[atom('O_4', 0)]) == {atom('C_4', 0), atom('C_b', 1)}
    assert graph.edges[atom('C_b', 1), atom('O_4', 0)] == {
        'bond_order': 1, 'linkage': 'beta-O-4',
    }
    for left, right, order in (
        ('C_1', 'C_2', 1), ('C_2', 'C_3', 2), ('C_3', 'C_4', 1),
        ('C_4', 'C_5', 1), ('C_5', 'C_6', 2), ('C_6', 'C_1', 1),
        ('C_4', 'O_4', 2), ('C_a', 'C_b', 1),
    ):
        assert graph.edges[atom(left, 1), atom(right, 1)]['bond_order'] == order
    assert graph.edges[atom('C_a', 2), atom('C_b', 2)]['bond_order'] == 1
    assert graph.edges[atom('C_a', 0), atom('C_b', 0)]['bond_order'] == 2
    assert not any(data['template_label'] in {'H_Oa', 'H_O4', 'O_b', 'H_Ob'}
                   for _, data in graph.nodes(data=True) if data['monomer_id'] == 1)
    assert atom('H_O4', 2) in graph
    assert nx.is_connected(graph)
    assert {data['monomer_id'] for _, data in graph.nodes(data=True)} == {0, 1, 2}
    for node, data in graph.nodes(data=True):
        assert data['index'] == node
        assert data['degree'] == graph.degree(node)
        valence = sum(graph.edges[node, other]['bond_order'] for other in graph[node])
        assert valence == {'C': 4, 'O': 2, 'H': 1}[data['element']]

    # Independent monolignol formulas: beta-ether formation adds net O;
    # spirodienone completion loses two H and adds no further oxygen.
    formulas = {'H': Counter(C=9, H=10, O=2), 'G': Counter(C=10, H=12, O=3),
                'S': Counter(C=11, H=14, O=4)}
    expected = sum((formulas[kind] for kind in monomers), Counter())
    expected['O'] += 1
    expected['H'] -= 2
    assert Counter(data['element'] for _, data in graph.nodes(data=True)) == expected


@pytest.mark.parametrize('precursor_linker,target', [('beta-O-4', 1), ('4-O-beta', 0)])
def test_spirodienone_reuses_precursor_oxygen_without_editing_inputs(tmp_path, precursor_linker, target):
    builder, chain = build(tmp_path, ['G', 'H'], [precursor_linker])
    incoming = builder.load_monomer('S', monomer_id=2, start_index=max(chain) + 1)
    originals = chain.copy(), incoming.copy()
    find = logic._find_attachment_atom
    oxygen = find(chain, 'O_a', target)
    removed = {find(chain, 'H_Oa', target), find(chain, 'H_O4', target)}
    result = logic.connect_spirodienone_beta_1(chain, incoming, target)

    assert nx.utils.graphs_equal(chain, originals[0])
    assert nx.utils.graphs_equal(incoming, originals[1])
    assert set(result) == (set(chain) | set(incoming)) - removed
    assert set(result[oxygen]) == {find(chain, 'C_a', target), find(incoming, 'C_a', 2)}
    assert sum(data.get('linkage') == SPIRO for _, _, data in result.edges(data=True)) == 3


@pytest.mark.parametrize('precursor_linker,target', [
    ('beta-O-4', 0), ('4-O-beta', 1), ('resinol-beta-beta', 1),
    ('phenylcoumaran-beta-5', 1),
])
def test_spirodienone_rejects_ineligible_existing_units(tmp_path, precursor_linker, target):
    builder, chain = build(tmp_path, ['G', 'H'], [precursor_linker])
    incoming = builder.load_monomer('G', monomer_id=2, start_index=max(chain) + 1)
    originals = chain.copy(), incoming.copy()
    with pytest.raises(ValueError, match='beta-O-4 precursor'):
        logic.connect_spirodienone_beta_1(chain, incoming, target)
    assert nx.utils.graphs_equal(chain, originals[0])
    assert nx.utils.graphs_equal(incoming, originals[1])


@pytest.mark.parametrize('linkers', [
    [SPIRO], ['4-O-beta', SPIRO], ['phenylcoumaran-beta-5', SPIRO],
    ['resinol-beta-beta', SPIRO], ['beta-O-4', SPIRO, SPIRO],
])
def test_spirodienone_recipe_requires_the_correct_precursor(tmp_path, linkers):
    with pytest.raises(ValueError, match="requires a preceding 'beta-O-4'"):
        build(tmp_path, ['G'] * (len(linkers) + 1), linkers)


@pytest.mark.parametrize('following', ['4-O-beta', 'resinol-beta-beta'])
def test_spirodienone_recipe_rejects_reusing_incoming_beta(tmp_path, following):
    with pytest.raises(ValueError, match='beta position of the shared monomer is already occupied'):
        build(tmp_path, ['G'] * 4, ['beta-O-4', SPIRO, following])


def test_legacy_beta1_name_is_rejected(tmp_path):
    with pytest.raises(ValueError, match='Unsupported linker'):
        build(tmp_path, ['G', 'H'], ['beta-1'])
