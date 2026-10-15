from types import SimpleNamespace

import numpy as np
import pytest

from MakroLyzer import dictionaries
from MakroLyzer.graphs import GraphManager
from MakroLyzer.dynamic_modules.RMSD import RMSDAnalyzer
from MakroLyzer.dynamic_modules.dynamicBase import OutputHandler
from MakroLyzer.dynamic_modules import dynamicAnalysisMain


def graph(coords, nodes=(0, 1), com=None):
    return SimpleNamespace(
        get_all_coordinates=lambda: (nodes, coords),
        get_com=lambda: np.mean(coords, axis=0) if com is None else com,
    )


@pytest.mark.parametrize('mode', ['collect', 'streaming', None])
def test_translation_and_output_with_reused_arrays(tmp_path, mode):
    path = tmp_path / 'rmsd.csv'
    handler = OutputHandler(path, mode) if mode else None
    analyzer = RMSDAnalyzer(handler)
    analyzer.initialize_output()
    coords = np.array([[-1., 0., 0.], [1., 0., 0.]])
    com = np.zeros(3)
    np.testing.assert_allclose(analyzer.run(graph(coords, com=com), 0), [0, 0])
    coords += [2, -3, 6]
    com[:] = [2, -3, 6]
    np.testing.assert_allclose(analyzer.run(graph(coords, com=com), 4), [7, 0])
    # Symmetric deformation survives COM correction, relative to frame zero.
    coords[0, 0] -= 1
    coords[1, 0] += 1
    np.testing.assert_allclose(analyzer.run(graph(coords, com=com), 8), [np.sqrt(50), 1])
    assert analyzer.frame_number == 3
    analyzer.finalize()
    analyzer.finalize_output()
    if handler:
        assert path.read_text().splitlines() == [
            'Frame, RMSD / Å, COM-corrected RMSD / Å',
            '0,0.000,0.000', '4,7.000,0.000', '8,7.071,1.000',
        ]


def test_mass_weighted_com_but_equal_atom_weighted_rmsd():
    molecular_graph = GraphManager()
    for node, element in enumerate(['C', 'H']):
        molecular_graph.add_node(node, element=element, x=0., y=0., z=0.)
    analyzer = RMSDAnalyzer()
    np.testing.assert_allclose(analyzer.compute(molecular_graph), [0, 0])
    molecular_graph.nodes[0]['x'] = 2.
    masses = dictionaries.dictMass()
    center = 2 * masses['C'] / (masses['C'] + masses['H'])
    np.testing.assert_allclose(analyzer.compute(molecular_graph),
                               [np.sqrt(2), np.sqrt(((2-center)**2 + center**2)/2)])


def test_rotation_is_not_removed():
    analyzer = RMSDAnalyzer()
    analyzer.compute(graph(np.array([[-1., 0, 0], [1., 0, 0]])))
    np.testing.assert_allclose(
        analyzer.compute(graph(np.array([[0., -1, 0], [0., 1, 0]]))),
        [np.sqrt(2), np.sqrt(2)],
    )


@pytest.mark.parametrize('coords,nodes', [
    (np.empty((0, 3)), ()),
    (np.zeros((2, 2)), (0, 1)),
    (np.zeros((1, 3)), (0, 1)),
    (np.array([[np.nan, 0, 0], [0, 0, 0]]), (0, 1)),
    (np.array([[np.inf, 0, 0], [0, 0, 0]]), (0, 1)),
])
def test_invalid_coordinates_do_not_initialize_reference(coords, nodes):
    analyzer = RMSDAnalyzer()
    with pytest.raises(ValueError, match='coordinates'):
        analyzer.run(graph(coords, nodes), 0)
    assert analyzer.refstruc is None
    assert analyzer.frame_number == 0


@pytest.mark.parametrize('com', [np.zeros(2), [np.nan, 0, 0], [np.inf, 0, 0]])
def test_invalid_com_does_not_initialize_reference(com):
    analyzer = RMSDAnalyzer()
    with pytest.raises(ValueError, match='center of mass'):
        analyzer.compute(graph(np.zeros((2, 3)), com=com))
    assert analyzer.refstruc is None


@pytest.mark.parametrize('nodes', [(1, 0), (0, 2), (0,)])
def test_changed_atom_selection_is_rejected(nodes):
    analyzer = RMSDAnalyzer()
    analyzer.run(graph(np.zeros((2, 3))), 0)
    with pytest.raises(ValueError, match='identities'):
        analyzer.run(graph(np.zeros((len(nodes), 3)), nodes), 1)
    assert analyzer.frame_number == 1
    np.testing.assert_allclose(analyzer.compute(graph(np.zeros((2, 3)))), [0, 0])


@pytest.mark.parametrize('static', [False, True])
def test_driver_unwraps_and_writes_both_rmsds(tmp_path, static):
    trajectory = tmp_path / 'trajectory.xyz'
    trajectory.write_text(''.join(f'1\nframe\nC {x} 0 0\n' for x in [8, 1, 4, 7, 0]))
    output = tmp_path / 'rmsd.csv'
    dynamicAnalysisMain.main({
        'xyzFile': str(trajectory), 'BoxSize': 10, 'nthStep': 2,
        'RMSD': True, 'RMSD_file': str(output), 'timestep': None,
        'staticTopology': static, 'subgraphSelection': ['C1'],
    })
    np.testing.assert_allclose(np.loadtxt(output, delimiter=',', skiprows=1),
                               [[0, 0, 0], [2, 6, 0], [4, 12, 0]])
