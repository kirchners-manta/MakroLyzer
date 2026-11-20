import csv
import sys

import pytest

from MakroLyzer import MakroLyzer
from MakroLyzer.graphs import GraphManager
from MakroLyzer.structure_modules.DistanceFromCOM import DistanceFromCOMAnalyzer
from MakroLyzer.structure_modules.structureBase import OutputHandler


def make_graph(*atoms):
    graph = GraphManager()
    for atom_id, element, x in atoms:
        graph.add_node(atom_id, element=element, x=x, y=0.0, z=0.0)
    return graph


@pytest.mark.parametrize('box, expected', [(None, 9.8), (10, 0.2), ([10, 12, 14], 0.2)])
def test_distance_uses_minimum_image_when_periodic(box, expected):
    graph = make_graph((1, 'C', 0.1), (2, 'O', 9.9))
    [(atom_ids, distance)] = DistanceFromCOMAnalyzer('C', 'O', box_size=box).compute(graph)
    assert atom_ids == {2}
    assert distance == pytest.approx(expected)


def test_solvent_com_unwraps_bonded_atoms():
    graph = make_graph((1, 'C', 0.5), (2, 'O', 9.8), (3, 'O', 0.2))
    graph.add_edge(2, 3)
    [(atom_ids, distance)] = DistanceFromCOMAnalyzer('C', 'O2', box_size=10).compute(graph)
    assert atom_ids == {2, 3}
    assert distance == pytest.approx(0.5)


def test_missing_particle_and_zero_mass():
    graph = make_graph((1, 'D', 0), (2, 'O', 2))
    with pytest.raises(ValueError, match='No particle matches'):
        DistanceFromCOMAnalyzer('C', 'O').compute(graph)
    with pytest.raises(ValueError, match='Total mass of the particle is zero'):
        DistanceFromCOMAnalyzer('D', 'O').compute(graph)


def test_csv_keeps_molecule_columns_across_frames(tmp_path):
    path = tmp_path / 'distances.csv'
    handler = OutputHandler(path, mode='streaming')
    analyzer = DistanceFromCOMAnalyzer('C', 'O', box_size=10, output_handler=handler)
    analyzer.run(make_graph((1, 'C', 0), (2, 'O', 2), (3, 'O', 3)), 0)
    assert path.read_text().splitlines() == ['Reference,C', 'O #1 (Å),O #2 (Å)', '0,2.000,3.000']
    analyzer.run(make_graph((1, 'C', 0), (3, 'O', 3), (2, 'O', 1), (4, 'O', 4)), 1)
    analyzer.run(make_graph((1, 'C', 0)), 2)
    analyzer.run(make_graph((1, 'C', 0), (3, 'O', 4)), 5)
    analyzer.finalize_output()

    with path.open() as output:
        assert list(csv.reader(output)) == [
            ['Reference', 'C'],
            ['O #1 (Å)', 'O #2 (Å)'],
            ['0', '2.000', '3.000'],
            ['1', '1.000', '3.000'],
            ['2', '', ''],
            ['5', '', '4.000'],
        ]


@pytest.mark.parametrize('filename', [None, 'custom.csv'])
def test_cli_distance_from_com(tmp_path, monkeypatch, filename):
    trajectory = tmp_path / 'trajectory.xyz'
    trajectory.write_text('2\nframe 0\nC 0.1 0 0\nO 8.1 0 0\n')
    output_dir = tmp_path / 'results'
    argv = [
        'MakroLyzer', '-xyz', str(trajectory), '--distance-from-com', 'C', 'O',
        '-bs', '10', '-sel', 'C', '--static-topology', '--output-dir', str(output_dir),
    ]
    if filename:
        argv.extend(['--distance-from-com-file', filename])
    monkeypatch.setattr(sys, 'argv', argv)

    MakroLyzer.main()

    with (output_dir / (filename or 'DistanceFromCOM.csv')).open() as output:
        assert list(csv.reader(output)) == [['Reference', 'C'], ['O #1 (Å)'], ['0', '2.000']]
