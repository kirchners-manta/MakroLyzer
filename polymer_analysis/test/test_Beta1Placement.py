"""Check spirodienone product geometry and continuation through user recipes."""

import json
from itertools import product

import numpy as np
import pytest

from MakroLyzer.builder_modules import liginConnectionLogic as logic
from MakroLyzer.builder_modules.ligninBuilder import LigninBuilder


SPIRO = 'spirodienone-beta-1'


# Check relaxed bond lengths for every H/G/S assignment, in addition to
# the chemical-connectivity/formula matrix in test_Beta1Topology.py.
@pytest.mark.parametrize('monomers,linkers', [
    (list(monomers), ['beta-O-4', SPIRO]) for monomers in product('HGS', repeat=3)
] + [
    (['G'] * 4, ['beta-O-4', 'beta-O-4', SPIRO]),
    (['G'] * 5, ['beta-O-4', SPIRO, 'beta-O-4', SPIRO]),
])
def test_spirodienone_has_reasonable_relaxed_bonds(tmp_path, monomers, linkers):
    recipe = tmp_path / 'trimer.json'
    recipe.write_text(json.dumps({'monomer': monomers, 'linker': linkers}))
    graph = LigninBuilder(recipe).run()
    find = logic._find_attachment_atom
    for middle, linker in enumerate(linkers):
        if linker != SPIRO:
            continue
        for left, left_id, right, right_id, expected in (
            ('C_1', middle, 'C_b', middle + 1, 1.54),
            ('C_a', middle, 'O_a', middle, 1.42),
            ('O_a', middle, 'C_a', middle + 1, 1.42),
            ('C_b', middle, 'O_4', middle - 1, 1.42),
        ):
            a, b = find(graph, left, left_id), find(graph, right, right_id)
            length = np.linalg.norm(logic._coordinates(graph, a) - logic._coordinates(graph, b))
            assert length == pytest.approx(expected, abs=0.15)

    # Check the entire product, including the retained ether partner. A
    # successful minimization alone can hide an accidentally strained bond.
    for left, right, data in graph.edges(data=True):
        elements = {graph.nodes[left]['element'], graph.nodes[right]['element']}
        length = np.linalg.norm(logic._coordinates(graph, left) - logic._coordinates(graph, right))
        if 'H' in elements:
            assert 0.85 < length < 1.2
        elif data['bond_order'] == 2:
            assert 1.1 < length < 1.5
        else:
            assert 1.25 < length < 1.75
    diagnostics = graph.graph['uff_optimization']
    assert np.isfinite([diagnostics['initial_energy'], diagnostics['final_energy']]).all()
    assert diagnostics['final_energy'] <= diagnostics['initial_energy']


@pytest.mark.parametrize('following', ['beta-O-4', 'phenylcoumaran-beta-5'])
def test_spirodienone_can_extend_at_incoming_units_free_ring_sites(tmp_path, following):
    recipe = tmp_path / 'extended.json'
    recipe.write_text(json.dumps({
        'monomer': ['G', 'S', 'G', 'H'], 'linker': ['beta-O-4', SPIRO, following],
    }))
    graph = LigninBuilder(recipe).run()
    assert {data['monomer_id'] for _, data in graph.nodes(data=True)} == {0, 1, 2, 3}
    find = logic._find_attachment_atom
    assert graph.has_edge(find(graph, 'C_1', 1), find(graph, 'C_b', 2))
    assert set(graph[find(graph, 'O_a', 1)]) == {find(graph, 'C_a', 1), find(graph, 'C_a', 2)}
    assert graph.edges[find(graph, 'C_b', 1), find(graph, 'O_4', 0)]['bond_order'] == 1


def test_spirodienone_recipe_runs_through_cli_and_writes_xyz(tmp_path, monkeypatch):
    from MakroLyzer import MakroLyzer

    recipe = tmp_path / 'cli.json'
    output = tmp_path / 'trimer.xyz'
    recipe.write_text(json.dumps({
        'monomer': ['G', 'H', 'S'], 'linker': ['beta-O-4', SPIRO],
    }))
    monkeypatch.setattr('sys.argv', [
        'MakroLyzer', '-buildLignin', str(recipe), '--lignin-file', str(output),
    ])
    MakroLyzer.main()
    rows = output.read_text().splitlines()
    # G/H/S trimer: C30H34O10, with one atom per XYZ row.
    assert int(rows[0]) == 74
    assert len(rows[2:]) == 74
