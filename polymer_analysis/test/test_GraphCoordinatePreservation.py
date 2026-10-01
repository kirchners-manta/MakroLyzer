"""Coordinate preservation must not change periodic topology or input graphs."""
import networkx as nx
import numpy as np
import pandas as pd
import pytest

from MakroLyzer.graphs import GraphManager


@pytest.mark.parametrize('box', [10., [10., 20., 30.]])
@pytest.mark.parametrize('direct', [False, True], ids=['constructor', 'create_graph'])
def test_preservation_keeps_all_axes_and_periodic_topology(box, direct):
    lengths = np.broadcast_to(box, (3,))
    # Different image offsets on every axis, including negative coordinates.
    coords = np.array([[.2, .2, .2], [1.6, .2, .2], [5., 5., 5.]])
    coords += np.array([[-2, 3, -1], [4, -2, 2], [-1, 1, 3]]) * lengths
    frame = pd.DataFrame(coords, columns=['x', 'y', 'z'])
    frame['atom'] = ['C', 'C', 'H']
    original = frame.copy(deep=True)
    if direct:
        graph = GraphManager()
        graph.create_graph(frame, boxSize=box, preserve_coords=True)
    else:
        graph = GraphManager(frame, boxSize=box, preserve_coords=True)
    np.testing.assert_array_equal(graph.get_all_coordinates()[1], coords)
    assert set(graph.edges) == {(0, 1)}
    assert [graph.nodes[i]['degree'] for i in graph] == [1, 1, 0]
    ordinary = GraphManager(frame, boxSize=box)
    assert set(ordinary.edges) == set(graph.edges)
    pd.testing.assert_frame_equal(frame, original)


@pytest.mark.parametrize('source_type', ['manager', 'networkx', 'selection'])
@pytest.mark.parametrize('preserve', [False, True])
def test_copy_preservation_and_default_shifting(source_type, preserve):
    source = GraphManager()
    # Nonconsecutive IDs check that selected graphs retain their identities.
    source.add_node(4, index=4, element='C', x=19.5, y=-20., z=30., degree=1)
    source.add_node(9, index=9, element='C', x=30.9, y=-20., z=30., degree=1)
    source.add_edge(4, 9, label='bond')
    if source_type == 'networkx':
        source = nx.Graph(source)
    elif source_type == 'selection':
        source.add_node(20, element='H', x=5., y=5., z=5.)
        source = source.subgraph([4, 9])
    original_nodes = {node: dict(attrs) for node, attrs in source.nodes(data=True)}
    result = GraphManager(source, boxSize=10, preserve_coords=preserve)
    assert list(result) == [4, 9]
    assert result.edges[4, 9]['label'] == 'bond'
    expected = [[19.5, -20., 30.], [30.9, -20., 30.]] if preserve else [[0., 0., 0.], [1.4, 0., 0.]]
    np.testing.assert_allclose(result.get_all_coordinates()[1], expected, atol=1e-12)
    assert dict(source.nodes(data=True)) == original_nodes
    # Coordinate edits on the copy must not affect the source graph either.
    result.nodes[4]['x'] = 999.
    assert source.nodes[4]['x'] == 19.5


@pytest.mark.parametrize('preserve', [False, True])
def test_without_box_preserves_coordinates_and_uses_euclidean_distances(preserve):
    frame = pd.DataFrame({'atom': ['C', 'C'], 'x': [19.5, 30.9],
                          'y': [-20., -20.], 'z': [30., 30.]})
    graph = GraphManager(frame, preserve_coords=preserve)
    np.testing.assert_array_equal(graph.get_all_coordinates()[1], frame[['x', 'y', 'z']])
    assert graph.number_of_edges() == 0


def test_preservation_does_not_override_custom_bond_cutoff():
    frame = pd.DataFrame({'atom': ['C', 'C'], 'x': [19.5, 30.9],
                          'y': [0., 0.], 'z': [0., 0.]})
    graph = GraphManager(frame, boxSize=10, vib_factor=.5, preserve_coords=True)
    assert graph.number_of_edges() == 0
    np.testing.assert_array_equal(graph.get_all_coordinates()[1], frame[['x', 'y', 'z']])


@pytest.mark.parametrize('box', [10., [10., 20., 30.]])
@pytest.mark.parametrize('preserve', [False, True])
def test_wrapped_and_unwrapped_frames_produce_same_topology(box, preserve):
    lengths = np.broadcast_to(box, (3,))
    # Follow a bonded carbon pair across the boundary over successive frames.
    # A distant hydrogen remains isolated, checking absent as well as present bonds.
    for displacement in [0., 1., 3., 7., 11.]:
        unwrapped = np.array([[9.5 + displacement, -1., 31.],
                              [10.9 + displacement, -1., 31.],
                              [5. + displacement, -5., 35.]])
        wrapped = unwrapped % lengths
        graphs = []
        for coords in (wrapped, unwrapped):
            frame = pd.DataFrame(coords, columns=['x', 'y', 'z'])
            frame['atom'] = ['C', 'C', 'H']
            graphs.append(GraphManager(frame, boxSize=box, preserve_coords=preserve))

        wrapped_graph, unwrapped_graph = graphs
        assert list(wrapped_graph.nodes) == list(unwrapped_graph.nodes)
        assert set(wrapped_graph.edges) == set(unwrapped_graph.edges) == {(0, 1)}
        for node in wrapped_graph:
            for attribute in ('index', 'element', 'degree'):
                assert wrapped_graph.nodes[node][attribute] == unwrapped_graph.nodes[node][attribute]

        wrapped_positions = wrapped_graph.get_all_coordinates()[1]
        unwrapped_positions = unwrapped_graph.get_all_coordinates()[1]
        if preserve:
            np.testing.assert_array_equal(wrapped_positions, wrapped)
            np.testing.assert_array_equal(unwrapped_positions, unwrapped)
        else:
            np.testing.assert_allclose(wrapped_positions, unwrapped_positions, atol=1e-12)
