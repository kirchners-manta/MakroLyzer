import pytest

from MakroLyzer.input_handling import checkInput
from MakroLyzer.input_handling import inputHandlingMain
from MakroLyzer.input_handling import readInput


def _base_args():
    return {
        'xyzFile': None,
        'lmpFile': None,
        'patternFile': None,
    }


class TestCheckInput:
    def test_xyz_file_valid(self, tmp_path):
        xyz_path = tmp_path / "traj.xyz"
        xyz_path.write_text("1\nC 0 0 0\n")

        args = _base_args()
        args['xyzFile'] = str(xyz_path)

        checkInput.checkInput(args)

    def test_xyz_file_missing(self):
        args = _base_args()
        args['xyzFile'] = "missing.xyz"

        with pytest.raises(checkInput.FileNotFoundError):
            checkInput.checkInput(args)

    def test_xyz_file_invalid_extension(self, tmp_path):
        bad_path = tmp_path / "traj.txt"
        bad_path.write_text("data")

        args = _base_args()
        args['xyzFile'] = str(bad_path)

        with pytest.raises(checkInput.InvalidFileFormatError):
            checkInput.checkInput(args)

    def test_xyz_file_empty(self, tmp_path):
        empty_path = tmp_path / "empty.xyz"
        empty_path.write_text("")

        args = _base_args()
        args['xyzFile'] = str(empty_path)

        with pytest.raises(checkInput.EmptyFileError):
            checkInput.checkInput(args)

    def test_pattern_file_valid(self, tmp_path):
        pattern_path = tmp_path / "patterns.txt"
        pattern_path.write_text("[[C,C]]\n")

        args = _base_args()
        args['patternFile'] = str(pattern_path)

        checkInput.checkInput(args)


class TestInputHandlingMain:
    def test_main_flags_analyzer_and_modifier(self, monkeypatch):
        args = {
            'xyzFile': "traj.xyz",
            'lmpFile': None,
            'nthStep': 1,
            'BoxSize': None,
            'hydrogenBonds': None,
            'hbondCube': None,
            'anisotropyFactor': None,
            'asphericityParameter': None,
            'radiusOfGyration': None,
            'NoSubgraphs': None,
            'endToEndDistance': None,
            'orderParameter': None,
            'formula': None,
            'Ramachandran': None,
            'patternFile': None,
            'saturation': None,
            'subgraph_coords': None,
        }

        monkeypatch.setattr(inputHandlingMain.readInput, 'readCommandLine', lambda: args)
        monkeypatch.setattr(inputHandlingMain.checkInput, 'checkInput', lambda a: None)

        analyzer, modifier, returned = inputHandlingMain.main(None)
        assert analyzer is False
        assert modifier is False
        assert returned == args

        args['anisotropyFactor'] = True
        args['patternFile'] = "patterns.txt"

        analyzer, modifier, returned = inputHandlingMain.main(None)
        assert analyzer is True
        assert modifier is True
        assert returned == args

    def test_main_exits_on_invalid_input(self, monkeypatch, capsys):
        args = _base_args()
        args['xyzFile'] = "missing.xyz"

        def _raise(_args):
            raise checkInput.FileNotFoundError("XYZ file 'missing.xyz' not found.")

        monkeypatch.setattr(inputHandlingMain.readInput, 'readCommandLine', lambda: args)
        monkeypatch.setattr(inputHandlingMain.checkInput, 'checkInput', _raise)

        with pytest.raises(SystemExit) as exc:
            inputHandlingMain.main(None)

        captured = capsys.readouterr()
        assert "XYZ file 'missing.xyz' not found." in captured.out
        assert exc.value.code == 1


class TestReadInputParsers:
    def test_ring_cycle_size_single_and_range(self):
        assert readInput.RingCycleSize("6") == 6
        assert readInput.RingCycleSize("[4,7]") == (4, 7)

    def test_ring_cycle_size_minimum_three(self):
        with pytest.raises(Exception, match="Minimum ring size is 3."):
            readInput.RingCycleSize("2")
        with pytest.raises(Exception, match="Minimum ring size is 3."):
            readInput.RingCycleSize("[2,7]")

    def test_hbond_cube_param(self):
        assert readInput.HBondCubeParam("100:20") == ((100.0, 100.0, 100.0), (20, 20, 20))
        assert readInput.HBondCubeParam("100,120,140:20,24,28") == (
            (100.0, 120.0, 140.0),
            (20, 24, 28),
        )
