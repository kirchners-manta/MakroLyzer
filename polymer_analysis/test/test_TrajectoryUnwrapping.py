"""Regression tests for continuous coordinates and periodic graph construction."""
import networkx as nx
import numpy as np
import pandas as pd
import pytest

from MakroLyzer.graphs import GraphManager
from MakroLyzer.input_handling.readLMP import readLMP
from MakroLyzer.input_handling.readXYZ import readXYZ


@pytest.fixture(params=['xyz', 'lmp'])
def trajectory(request, tmp_path):
    """Two atoms crossing repeatedly in opposite directions, on all axes."""
    box = np.array([10., 20., 30.])
    expected = np.array([[.9, .1], [1.1, -.1], [1.3, -.3],
                         [1.7, -.7], [2.1, -1.1]])[:, :, None] * box
    wrapped = expected % box
    path = tmp_path / request.param
    blocks = []
    for step, positions in enumerate(wrapped):
        if request.param == 'xyz':
            blocks.append('2\nframe\n' + ''.join(
                f'{element} {x} {y} {z}\n'
                for element, (x, y, z) in zip(['C', 'H'], positions)))
        else:
            blocks.append('ITEM: TIMESTEP\n' + str(step) +
                          '\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS pp pp pp\n'
                          '0 10\n0 20\n0 30\nITEM: ATOMS id element x y z\n' +
                          ''.join(f'{i + 1} {element} {x} {y} {z}\n'
                                  for i, (element, (x, y, z)) in
                                  enumerate(zip(['C', 'H'], positions))))
    path.write_text(''.join(blocks))
    reader = readXYZ if request.param == 'xyz' else readLMP
    return reader, path, box, expected, wrapped


def coordinates(frames):
    return np.array([frame[['x', 'y', 'z']].to_numpy() for frame in frames])


def test_repeated_crossings(trajectory):
    reader, path, box, expected, _ = trajectory
    np.testing.assert_allclose(coordinates(reader(path, unwrap=True, box_size=box)), expected)


def test_disabled_unwrapping(trajectory):
    reader, path, _, _, wrapped = trajectory
    np.testing.assert_allclose(coordinates(reader(path)), wrapped)


def test_caller_edits_do_not_change_unwrapping_state(trajectory):
    reader, path, box, expected, _ = trajectory
    frames = reader(path, unwrap=True, box_size=box)
    first = next(frames)
    first.loc[:, ['x', 'y', 'z']] = 999.
    np.testing.assert_allclose(coordinates(frames), expected[1:])


@pytest.mark.parametrize('box', [None, 0, -1, [10, 20], [10, np.nan, 30]])
def test_invalid_box(trajectory, box):
    reader, path, _, _, _ = trajectory
    with pytest.raises(ValueError, match='box_size'):
        list(reader(path, unwrap=True, box_size=box))


def write_dump(path, frames, columns='id element x y z'):
    path.write_text(''.join(
        f'ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n{len(rows)}\n'
        'ITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n'
        f'ITEM: ATOMS {columns}\n' + '\n'.join(rows) + '\n'
        for step, rows in enumerate(frames)))


def test_lmp_matches_reordered_atom_ids(tmp_path):
    path = tmp_path / 'reordered.lmp'
    write_dump(path, [['1 C 9 0 0', '2 H 1 0 0'],
                      ['2 H 9 0 0', '1 C 1 0 0']])
    frames = list(readLMP(path, unwrap=True, box_size=10))
    assert frames[1]['id'].tolist() == ['1', '2']
    np.testing.assert_allclose(frames[1].x, [11, -1])


@pytest.mark.parametrize('unwrap', [False, True])
def test_lmp_native_unwrapped_takes_precedence_without_box(tmp_path, unwrap):
    path = tmp_path / 'unwrapped.lmp'
    # The large real displacement must not receive minimum-image correction.
    write_dump(path, [['1 C 9 0 0 19 0 0'], ['1 C 1 0 0 41 0 0']],
               'id element x y z xu yu zu')
    assert [f.x.iloc[0] for f in readLMP(path, unwrap=unwrap)] == [19, 41]


def test_lmp_image_counters_are_not_coordinates(tmp_path):
    path = tmp_path / 'images.lmp'
    write_dump(path, [['1 C 0 0 0']], 'id element ix iy iz')
    with pytest.raises(ValueError, match='coordinate columns'):
        list(readLMP(path))


@pytest.mark.parametrize('rows, message', [
    (['1 C 1 0 0'], 'Atom count'),
    (['1 C 1 0 0', '3 H 9 0 0'], 'Atom IDs'),
    (['1 C 1 0 0', '1 H 9 0 0'], 'Duplicate'),
])
def test_lmp_rejects_changed_identity(tmp_path, rows, message):
    path = tmp_path / 'identity.lmp'
    write_dump(path, [['1 C 9 0 0', '2 H 1 0 0'], rows])
    with pytest.raises(ValueError, match=message):
        list(readLMP(path, unwrap=True, box_size=10))


@pytest.mark.parametrize('preserve', [False, True])
def test_graph_preserves_coordinates_and_periodic_bonds(preserve):
    # Atoms are in different periodic images but are bonded by minimum image.
    frame = pd.DataFrame({'atom': ['C', 'C'], 'x': [19.5, 30.9],
                          'y': [0., 0.], 'z': [0., 0.]})
    graph = GraphManager(frame, boxSize=10, preserve_coords=preserve)
    assert graph.has_edge(0, 1)
    expected = [19.5, 30.9] if preserve else [0, 1.4]
    np.testing.assert_allclose(graph.get_all_coordinates()[1][:, 0], expected, atol=1e-12)
    for source in (graph, nx.Graph(graph), graph.subgraph([0, 1])):
        copy = GraphManager(source, boxSize=10, preserve_coords=preserve)
        assert copy.has_edge(0, 1)
        np.testing.assert_allclose(copy.get_all_coordinates()[1][:, 0], expected, atol=1e-12)
    np.testing.assert_allclose(frame.x, [19.5, 30.9])
