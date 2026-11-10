"""Verify optimization limits from the CLI reach each builder step."""

import json
from unittest.mock import Mock

import networkx as nx
import pytest

from MakroLyzer.builder_modules import builder_registry, globalHelpers
from MakroLyzer.input_handling import readInput


@pytest.mark.parametrize('limit', [None, '25'])
def test_iteration_limit_reaches_each_optimization(monkeypatch, tmp_path, limit):
    recipe = tmp_path / 'recipe.json'
    recipe.write_text(json.dumps({
        'monomer': ['G', 'H', 'S'], 'linker': ['beta-O-4', 'beta-O-4'],
    }))
    argv = ['MakroLyzer', '-buildLignin', str(recipe)]
    if limit is not None:
        argv += ['--builder-max-iterations', limit]
    monkeypatch.setattr('sys.argv', argv)
    builder = builder_registry.create_lignin(readInput.readCommandLine())
    monkeypatch.setattr(builder, '_connect_monomers',
                        lambda chain, linker, monomer, **kw: nx.compose(chain, monomer))
    optimize = Mock(side_effect=lambda graph, **kw: graph)
    monkeypatch.setattr(globalHelpers, 'optimize_structure', optimize)
    builder.build_structure()
    assert optimize.call_count == 2
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
