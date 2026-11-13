"""Verify optimization limits from the CLI reach each builder step."""

import json
from unittest.mock import Mock

import networkx as nx
import pytest

from MakroLyzer.builder_modules import builder_registry, globalHelpers
from MakroLyzer.input_handling import readInput
from MakroLyzer.builder_modules.ligninBuilder import LigninBuilder


@pytest.mark.parametrize('monomers,linkers,blocked', [
    (['G', 'S'], ['beta-5'], True),
    (['H', 'G', 'S'], ['4-O-beta', 'beta-5'], True),
    (['S', 'G'], ['beta-5'], False),
    (['G', 'H'], ['beta-5'], True),
    (['H', 'G'], ['beta-5'], False),
    (['G', 'S'], ['4-O-beta'], False),
])
def test_monomer_linkage_rules_are_directional(monkeypatch, monomers, linkers, blocked):
    # Restrict the existing monomer/linker pair regardless of incoming type.
    monkeypatch.setattr(LigninBuilder, '_forbidden_monomer_linkages', staticmethod(
        lambda: {('G', 'beta-5'): 'Test restriction'}
    ))
    builder = LigninBuilder.__new__(LigninBuilder)
    builder.recipe = {'monomer': monomers, 'linker': linkers}
    if blocked:
        with pytest.raises(ValueError, match='Forbidden monomer linkage.*Test restriction'):
            builder.validate_recipe()
    else:
        builder.validate_recipe()


@pytest.mark.parametrize('monomers', [
    ['G', 'H', 'S'], ['S', 'G', 'H'], ['H', 'S', 'G'],
])
@pytest.mark.parametrize('linkers,site', [
    (['beta-O-4', '4-O-beta'], 'beta'),
    (['4-O-beta', 'beta-O-4'], 'O-4'),
])
def test_consecutive_linkers_cannot_reuse_middle_site(monomers, linkers, site):
    builder = LigninBuilder.__new__(LigninBuilder)
    builder.recipe = {'monomer': monomers, 'linker': linkers}
    with pytest.raises(ValueError, match=f'{site} position'):
        builder.validate_recipe()


@pytest.mark.parametrize('linkers', [
    [], ['beta-O-4'], ['4-O-beta'],
    ['beta-O-4', 'beta-O-4'], ['4-O-beta', '4-O-beta'],
])
def test_distinct_middle_sites_are_allowed(linkers):
    builder = LigninBuilder.__new__(LigninBuilder)
    builder.recipe = {'monomer': ['G'] * (len(linkers) + 1), 'linker': linkers}
    builder.validate_recipe()


@pytest.mark.parametrize('limit', [None, '25'])
def test_iteration_limit_reaches_each_optimization(monkeypatch, tmp_path, limit):
    recipe = tmp_path / 'recipe.json'
    recipe.write_text(json.dumps({
        'monomer': ['G', 'H', 'S'], 'linker': ['4-O-beta', '4-O-beta'],
    }))
    argv = ['MakroLyzer', '-buildLignin', str(recipe)]
    if limit is not None:
        argv += ['--builder-max-iterations', limit]
    monkeypatch.setattr('sys.argv', argv)
    builder = builder_registry.create_lignin(readInput.readCommandLine())
    monkeypatch.setattr(builder, '_connect_monomers',
                        lambda chain, linker, monomer, **kw: nx.compose(chain, monomer))
    def optimized_copy(graph, **kwargs):
        result = graph.copy()
        result.graph['uff_optimization'] = {'final_energy': 10.0}
        return result

    optimize = Mock(side_effect=optimized_copy)
    monkeypatch.setattr(globalHelpers, 'optimize_structure', optimize)
    builder.build_structure()
    assert optimize.call_count == 3
    assert all(call.kwargs == {'max_iterations': int(limit) if limit else 1500}
               for call in optimize.call_args_list)


@pytest.mark.parametrize('limit', ['0', '-1', '1.5', 'abc'])
def test_invalid_iteration_limit_rejected(monkeypatch, limit):
    monkeypatch.setattr('sys.argv', [
        'MakroLyzer', '-buildLignin', 'recipe.json', '--builder-max-iterations', limit,
    ])
    with pytest.raises(SystemExit) as exc:
        readInput.readCommandLine()
    assert exc.value.code == 2


@pytest.mark.parametrize('second_energy,keep_second', [(5.0, True), (15.0, False), (10.0, False)])
def test_final_relaxation_keeps_lower_energy_result(monkeypatch, second_energy, keep_second):
    builder = LigninBuilder.__new__(LigninBuilder)
    builder.recipe = {'monomer': ['G', 'H'], 'linker': ['beta-5']}
    builder.max_iterations = 25
    monkeypatch.setattr(builder, '_connect_monomers',
                        lambda chain, linker, monomer, **kw: nx.compose(chain, monomer))
    results = []

    def optimize(graph, **kwargs):
        result = graph.copy()
        energy = 10.0 if not results else second_energy
        result.graph['uff_optimization'] = {'final_energy': energy}
        results.append(result)
        return result

    monkeypatch.setattr(globalHelpers, 'optimize_structure', optimize)
    result = builder.build_structure()
    assert len(results) == 2
    assert result is results[1 if keep_second else 0]
    assert results[0].graph['uff_optimization']['final_energy'] == 10.0


def test_single_monomer_receives_one_optimization(monkeypatch):
    builder = LigninBuilder.__new__(LigninBuilder)
    builder.recipe = {'monomer': ['G'], 'linker': []}
    builder.max_iterations = 25
    optimize = Mock(side_effect=lambda graph, **kw: graph)
    monkeypatch.setattr(globalHelpers, 'optimize_structure', optimize)
    builder.build_structure()
    assert optimize.call_count == 1
    assert optimize.call_args.kwargs == {'max_iterations': 25}
