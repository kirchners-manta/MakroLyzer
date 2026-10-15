from types import SimpleNamespace

import numpy as np
import pytest

from MakroLyzer.dynamic_modules.MSD import MSDAnalyzer
from MakroLyzer.dynamic_modules.dynamicBase import OutputHandler
from MakroLyzer.dynamic_modules import dynamicAnalysisMain


def graph(coords, nodes=(0, 1)):
    return SimpleNamespace(get_all_coordinates=lambda: (nodes, coords))


@pytest.mark.parametrize('mode', ['collect', 'streaming', None])
def test_average_over_atoms_and_origins_with_reused_coordinates(tmp_path, mode):
    path = tmp_path / 'msd.csv'
    handler = OutputHandler(path, mode) if mode else None
    analyzer = MSDAnalyzer(1.1, .5, handler)
    analyzer.initialize_output()
    coords = np.zeros((2, 3))
    for i, x in enumerate([0, 1, 3, 6]):
        coords[:, 0] = [x, 2*x]
        analyzer.run(graph(coords), i)
    # Per-origin values at lag 1 are 2.5, 10, 22.5; lag 2: 22.5, 62.5.
    results = analyzer.finalize()
    np.testing.assert_allclose(results['lag_time'], [.5, 1.])
    np.testing.assert_allclose(results['msd'], [35/3, 42.5])
    np.testing.assert_array_equal(results['n_origins'], [3, 2])
    assert len(analyzer.structures) == 2
    assert analyzer.frame_number == 4
    if handler:
        assert handler.accumulated_rows == []
    analyzer.finalize_output()
    if handler:
        rows = path.read_text().splitlines()
        assert rows[0] == 'Correlation Time, MSD / Å²'
        assert len(rows) == 3
        np.testing.assert_allclose(np.loadtxt(path, delimiter=',', skiprows=1)[:, 1], results['msd'])


@pytest.mark.parametrize('n_frames', [0, 1])
def test_short_trajectory_has_no_sampled_lags(n_frames):
    analyzer = MSDAnalyzer(2, 1)
    for i in range(n_frames):
        analyzer.run(graph(np.zeros((2, 3))), i)
    assert analyzer.finalize()['msd'].size == 0


def test_changed_selection_is_rejected():
    analyzer = MSDAnalyzer(2, 1)
    analyzer.run(graph(np.zeros((2, 3))), 0)
    with pytest.raises(ValueError, match='identities'):
        analyzer.run(graph(np.zeros((2, 3)), nodes=(1, 0)), 1)


@pytest.mark.parametrize('depth, timestep', [(0, 1), (1, 0), (-1, 1), (1, np.nan)])
def test_invalid_timing(depth, timestep):
    with pytest.raises(ValueError):
        MSDAnalyzer(depth, timestep)


@pytest.mark.parametrize('static', [False, True])
def test_driver_unwraps_before_stride_and_finalizes(tmp_path, static):
    path = tmp_path / 'trajectory.xyz'
    # True motion is +3 per saved frame. Applying no-jump after stride would
    # mistake +6 for -4 in this box.
    path.write_text(''.join(f'1\nframe\nC {x} 0 0\n' for x in [8, 1, 4, 7, 0]))
    output = tmp_path / 'msd.csv'
    dynamicAnalysisMain.main({
        'xyzFile': str(path), 'BoxSize': 10, 'nthStep': 2,
        'timestep': .5, 'MSD': 2, 'MSD_file': str(output),
        'staticTopology': static, 'subgraphSelection': ['C1'],
    })
    np.testing.assert_allclose(np.loadtxt(output, delimiter=',', skiprows=1),
                               [[1, 36], [2, 144]])


def analyze_positions(positions, depth, timestep=1.):
    """Run a trajectory with a fixed identity for each atom."""
    analyzer = MSDAnalyzer(depth, timestep)
    nodes = tuple(range(positions.shape[1]))
    for frame_idx, coords in enumerate(positions):
        analyzer.run(graph(coords, nodes), frame_idx)
    return analyzer.finalize()


def test_stationary_atoms_have_zero_msd_at_every_lag():
    positions = np.repeat([[[2., -3., 4.], [-8., 7., 1.]]], 6, axis=0)
    result = analyze_positions(positions, depth=10)
    np.testing.assert_array_equal(result['msd'], np.zeros(5))
    np.testing.assert_array_equal(result['n_origins'], [5, 4, 3, 2, 1])
    np.testing.assert_array_equal(result['lag_time'], [1, 2, 3, 4, 5])


def test_constant_velocity_gives_quadratic_msd_in_physical_time():
    # |v|² = 49. All three axes contribute, without dividing by dimension.
    velocity = np.array([2., -3., 6.])
    times = np.arange(8) * .25
    positions = np.array([[[5., 1., -2.], [-3., 4., 8.]]]) + times[:, None, None] * velocity
    result = analyze_positions(positions, depth=1., timestep=.25)
    np.testing.assert_allclose(result['msd'], 49 * result['lag_time']**2)
    np.testing.assert_array_equal(result['n_origins'], [7, 6, 5, 4])


def test_opposite_atom_motion_does_not_cancel_like_center_of_mass_motion():
    positions = np.array([[[0., 0., 0.], [0., 0., 0.]],
                          [[1., 2., 2.], [-1., -2., -2.]]])
    result = analyze_positions(positions, depth=1)
    # Each atom moves a distance of 3, although the center of mass is fixed.
    np.testing.assert_allclose(result['msd'], [9.])


@pytest.mark.parametrize('depth, timestep, expected_lags', [
    (.09, .1, []),          # No positive sampled lag fits the window.
    (.29, .1, [.1, .2]),    # Do not round the window up.
    (.3, .1, [.1, .2, .3]), # Include the boundary despite float roundoff.
    (10., .5, [.5, 1., 1.5]), # Do not invent samples beyond the trajectory.
])
def test_time_window_includes_only_available_sampled_lags(depth, timestep, expected_lags):
    positions = np.zeros((4, 1, 3))
    positions[:, 0, 0] = np.arange(4)
    result = analyze_positions(positions, depth, timestep)
    np.testing.assert_allclose(result['lag_time'], expected_lags)
    np.testing.assert_allclose(result['msd'], np.arange(1, len(expected_lags) + 1)**2)
    np.testing.assert_array_equal(result['n_origins'], 4 - np.arange(1, len(expected_lags) + 1))


@pytest.mark.parametrize('n_atoms', [1, 5])
@pytest.mark.parametrize('max_lag', [1, 4, 11])
def test_matches_independent_all_origin_reference(n_atoms, max_lag):
    # Deterministic irregular motion exercises averaging after history eviction.
    rng = np.random.default_rng(182)
    positions = np.cumsum(rng.normal(size=(12, n_atoms, 3)), axis=0)
    result = analyze_positions(positions, depth=max_lag * .2, timestep=.2)
    expected = []
    for lag in range(1, max_lag + 1):
        squared_distances = []
        for start in range(len(positions) - lag):
            for atom in range(n_atoms):
                squared_distances.append(sum(
                    (positions[start + lag, atom, axis] - positions[start, atom, axis])**2
                    for axis in range(3)
                ))
        expected.append(sum(squared_distances) / len(squared_distances))
    np.testing.assert_allclose(result['msd'], expected)
    np.testing.assert_allclose(result['lag_time'], np.arange(1, max_lag + 1) * .2)
    np.testing.assert_array_equal(result['n_origins'], 12 - np.arange(1, max_lag + 1))


def test_msd_is_invariant_to_fixed_atom_offsets_and_rotation():
    rng = np.random.default_rng(42)
    positions = rng.normal(size=(7, 3, 3))
    angle = .7
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0.],
                         [np.sin(angle), np.cos(angle), 0.],
                         [0., 0., 1.]])
    # Each atom may start in a different periodic image; constant offsets cancel.
    offsets = np.array([[20., 0., -10.], [0., 30., 10.], [-20., 10., 0.]])
    original = analyze_positions(positions, depth=3)
    transformed = analyze_positions(positions @ rotation.T + offsets, depth=3)
    np.testing.assert_allclose(transformed['msd'], original['msd'])


def test_finalizing_mid_trajectory_does_not_reset_or_duplicate_samples():
    analyzer = MSDAnalyzer(2, 1)
    for i, x in enumerate([0., 1.]):
        analyzer.run(graph(np.array([[x, 0., 0.]]), nodes=(0,)), i)
    first = analyzer.finalize()
    np.testing.assert_allclose(first['msd'], [1.])
    np.testing.assert_allclose(analyzer.finalize()['msd'], [1.])
    analyzer.run(graph(np.array([[3., 0., 0.]]), nodes=(0,)), 2)
    result = analyzer.finalize()
    np.testing.assert_allclose(result['msd'], [2.5, 9.])
    np.testing.assert_array_equal(result['n_origins'], [2, 1])
    np.testing.assert_allclose(first['msd'], [1.])
