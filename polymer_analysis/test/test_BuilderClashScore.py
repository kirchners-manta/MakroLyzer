"""Compare the spatial clash search with an exhaustive reference."""

import networkx as nx
import numpy as np
import pytest
from unittest.mock import patch

from MakroLyzer import dictionaries
from MakroLyzer.builder_modules import liginConnectionLogic as logic
from MakroLyzer.builder_modules.ligninBuilder import LigninBuilder


def make_graph(positions, prefix):
    graph = nx.Graph()
    for i, position in enumerate(positions):
        graph.add_node(
            f'{prefix}{i}', element=('C', 'H', 'O')[i % 3],
            x=position[0], y=position[1], z=position[2],
        )
    return graph


@pytest.mark.parametrize('scale', [0.5, 0.75, 1.0])
def test_kdtree_matches_exhaustive_pairs(scale):
    rng = np.random.default_rng(2)
    fixed = make_graph(rng.normal(size=(100, 3)) * 4, 'a')
    moving = make_graph(rng.normal(size=(20, 3)) * 4, 'b')
    # Coincident atoms exercise an actual exclusion, not a distant unused pair.
    moving.nodes['b0'].update({key: fixed.nodes['a0'][key] for key in ('x', 'y', 'z')})
    exclusions = {('a0', 'b0'), ('a4', 'b5')}
    radii = dictionaries.dictVdW()
    expected = 0.0
    for left, a in fixed.nodes(data=True):
        for right, b in moving.nodes(data=True):
            if (left, right) in exclusions:
                continue
            distance = np.linalg.norm(
                np.array([a[k] for k in ('x', 'y', 'z')])
                - np.array([b[k] for k in ('x', 'y', 'z')])
            )
            expected += max(scale * (radii[a['element']] + radii[b['element']]) - distance, 0) ** 2
    assert logic.clash_score(fixed, moving, exclusions, scale) == pytest.approx(expected)


def test_empty_and_distant_graphs_have_zero_score():
    fixed = make_graph([(0, 0, 0)], 'a')
    assert logic.clash_score(fixed, nx.Graph()) == 0
    assert logic.clash_score(fixed, make_graph([(100, 0, 0)], 'b')) == 0


def test_zero_score_stops_linker_search():
    builder = LigninBuilder.__new__(LigninBuilder)
    chain = builder.load_monomer('G')
    incoming = builder.load_monomer('H', monomer_id=1, start_index=25)
    with patch.object(logic, 'clash_score', return_value=0) as score:
        result = logic.connect_beta_O_4(chain, incoming, 0)
    assert score.call_count == 1
    assert len(result) == 47
    assert score.call_args.kwargs['excluded_pairs']
