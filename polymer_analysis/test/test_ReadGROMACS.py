import sys
from unittest.mock import Mock

import MDAnalysis as mda
import numpy as np
import pytest

from MakroLyzer.input_handling import checkInput, readInput
from MakroLyzer.input_handling.estimateFrames import EstimateFrames
from MakroLyzer.input_handling.readGROMACS import readGROMACS


@pytest.fixture
def gro(tmp_path):
    path = tmp_path / 'water.gro'
    lines = ['Water', '3']
    for index, (name, x) in enumerate([('OW', 0.1), ('HW1', 0.196), ('HW2', 0.068)], 1):
        lines.append(f'{1:5d}{"SOL":<5}{name:>5}{index:5d}{x:8.3f}{0.2:8.3f}{0.3:8.3f}')
    path.write_text('\n'.join(lines + ['   2.00000   2.00000   2.00000']) + '\n')
    return path


def test_gro_elements_units_and_indices(gro):
    assert EstimateFrames.estimateFramesGROMACS(gro) == 1
    with pytest.warns(UserWarning, match='guessing elements'):
        frames = list(readGROMACS(gro))
    assert len(frames) == 1
    assert frames[0].columns.tolist() == ['atom', 'x', 'y', 'z', 'index']
    assert frames[0]['atom'].tolist() == ['O', 'H', 'H']
    assert frames[0]['index'].tolist() == [0, 1, 2]
    np.testing.assert_allclose(frames[0].loc[0, ['x', 'y', 'z']].to_numpy(dtype=float), [1, 2, 3])


@pytest.mark.parametrize('extension', ['trr', 'xtc', 'pdb'])
def test_multiple_frames_and_independent_coordinates(gro, tmp_path, extension):
    universe = mda.Universe(str(gro), dt=1.0)
    universe.guess_TopologyAttrs(to_guess=['elements'])
    if extension == 'pdb':
        # GRO has no PDB metadata; supply explicit values for the test writer.
        for attribute, values in {
            'altLocs': [''] * 3,
            'icodes': [''],
            'chainIDs': ['A'] * 3,
            'occupancies': [1.0] * 3,
            'tempfactors': [0.0] * 3,
            'record_types': ['ATOM'] * 3,
            'formalcharges': [0] * 3,
        }.items():
            universe.add_TopologyAttr(attribute, values)
    path = tmp_path / f'trajectory.{extension}'
    original = universe.atoms.positions.copy()
    with mda.Writer(str(path), n_atoms=3, multiframe=True) as writer:
        writer.write(universe.atoms)
        universe.atoms.positions += 1.0
        universe.trajectory.ts.time = 1.0
        writer.write(universe.atoms)
    universe.trajectory.close()
    assert EstimateFrames.estimateFramesGROMACS(path) == 2
    if extension == 'pdb':
        frames = list(readGROMACS(path))
    else:
        with pytest.warns(UserWarning, match='guessing elements'):
            frames = list(readGROMACS(gro, path))
    assert len(frames) == 2
    assert frames[0]['atom'].tolist() == ['O', 'H', 'H']
    np.testing.assert_allclose(frames[0][['x', 'y', 'z']], original, atol=1e-3)
    np.testing.assert_allclose(frames[1][['x', 'y', 'z']], original + 1, atol=1e-3)


@pytest.mark.parametrize('flag', ['gro', 'pdb', 'trr', 'xtc'])
def test_cli(flag, monkeypatch):
    command = ['MakroLyzer', f'-{flag}', f'input.{flag}']
    if flag in ('trr', 'xtc'):
        command += ['-tpr', 'topol.tpr']
    monkeypatch.setattr(sys, 'argv', command)
    args = readInput.readCommandLine()
    assert args[f'{flag}File'] == f'input.{flag}'
    assert args['topologyFile'] == ('topol.tpr' if flag in ('trr', 'xtc') else None)


def test_cli_rejects_multiple_inputs(monkeypatch):
    monkeypatch.setattr(sys, 'argv', ['MakroLyzer', '-xyz', 'a.xyz', '-gro', 'a.gro'])
    with pytest.raises(SystemExit):
        readInput.readCommandLine()


def test_input_validation(gro, tmp_path):
    args = dict(xyzFile=None, lmpFile=None, patternFile=None, groFile=str(gro))
    checkInput.checkInput(args)
    args.update(groFile=None, trrFile=str(tmp_path / 'traj.trr'))
    with pytest.raises(ValueError, match='requires -tpr'):
        checkInput.checkInput(args)
    args['topologyFile'] = str(gro)
    with pytest.raises(checkInput.FileNotFoundError):
        checkInput.checkInput(args)
    path = tmp_path / 'traj.trr'
    path.touch()
    with pytest.raises(checkInput.EmptyFileError):
        checkInput.checkInput(args)
    path.write_bytes(b'placeholder')
    args['topologyFile'] = str(path)
    with pytest.raises(checkInput.InvalidFileFormatError):
        checkInput.checkInput(args)
    args.update(trrFile=None, topologyFile=str(gro), groFile=str(gro))
    with pytest.raises(ValueError, match='requires a TRR'):
        checkInput.checkInput(args)


@pytest.mark.parametrize('pipeline', ['analysis', 'modification'])
def test_gromacs_reaches_processing_pipeline(gro, tmp_path, monkeypatch, pipeline):
    from MakroLyzer.structure_modules import structureAnalysisMain
    from MakroLyzer.modify_modules import structureModificationMain

    universe = mda.Universe(str(gro), dt=1.0)
    path = tmp_path / 'trajectory.trr'
    with mda.Writer(str(path), n_atoms=3) as writer:
        for frame in range(3):
            universe.trajectory.ts.time = float(frame)
            writer.write(universe.atoms)
    universe.trajectory.close()
    module, registry = (
        (structureAnalysisMain, 'ANALYZERS_REGISTRATION') if pipeline == 'analysis'
        else (structureModificationMain, 'MODIFIERS_REGISTRATION')
    )
    processor = Mock(requires_full_graph=False)
    monkeypatch.setattr(module, registry, {'test': lambda args, **kwargs: processor})
    progress = Mock(side_effect=lambda frames, **kwargs: frames)
    monkeypatch.setattr(module, 'tqdm', progress)
    with pytest.warns(UserWarning, match='guessing elements'):
        module.main(dict(trrFile=str(path), topologyFile=str(gro), nthStep=2, orderParameter=None))
    assert progress.call_args.kwargs['total'] == 3
    assert [call.args[1] for call in processor.run.call_args_list] == [0, 2]
    graph = processor.run.call_args.args[0]
    assert len(graph) == 3
    assert graph.number_of_edges() == 2
