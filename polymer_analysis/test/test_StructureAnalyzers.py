import pytest
import tempfile
from pathlib import Path
import numpy as np

from unittest.mock import Mock, patch

from src.MakroLyzer.input_handling import readXYZ
from src.MakroLyzer.structure_modules import graphs
from src.MakroLyzer.errorOutputs.ErrorOutputs import ErrorOutputs

from src.MakroLyzer.structure_modules.structureBase import OutputHandler, StructureAnalyzer
from src.MakroLyzer.structure_modules.Hbonds import HBondsAnalyzer
from src.MakroLyzer.structure_modules.Anisotropy import AnisotropyAnalyzer
from src.MakroLyzer.structure_modules.Asphericity import AsphericityAnalyzer
from src.MakroLyzer.structure_modules.RadiusOfGyration import RadiusOfGyrationAnalyzer
from src.MakroLyzer.structure_modules.MoleculeCount import MoleculeCountAnalyzer
from src.MakroLyzer.structure_modules.EndToEndDistance import EndToEndDistanceAnalyzer
from src.MakroLyzer.structure_modules.OrderParameter import OrderParameterAnalyzer
from src.MakroLyzer.structure_modules.Ramachandran import RamachandranAnalyzer
from src.MakroLyzer.structure_modules.Dihedrals import DihedralsAnalyzer
from src.MakroLyzer.structure_modules.ChemicalFormula import ChemicalFormulaAnalyzer

# OutputHandler Tests # ------------------------------------------------------------------
class TestOutputHandler:
    """Test the OutputHandler class."""
    
    @pytest.fixture
    def temp_file(self):
        """Fixture that creates a temporary file path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir) / "output.csv"
            
    def test_initialization(self, temp_file):
        """Test initialization of OutputHandler."""
        handler = OutputHandler(temp_file, mode='collect')
        assert handler.file_path == temp_file
        assert handler.mode == 'collect'
        assert handler.accumulated_rows == []
        assert handler._initialized is False
    
    def test_write_csv(self, temp_file):
        """Test writing to a CSV file."""
        handler = OutputHandler(temp_file, mode='collect')
        header = "Column1, Column2"
        rows = ["1, 2", "3, 4"]
        handler.write_csv(header, rows)
        
        assert temp_file.exists()
        content = temp_file.read_text()
        assert content == "Column1, Column2\n1, 2\n3, 4\n"
        
    def test_initialize_file(self, temp_file):
        """Test initializing a file in streaming mode."""
        handler = OutputHandler(temp_file, mode='streaming')
        header = "Column1, Column2, Column3"
        handler.initialize_file(header)
        
        assert temp_file.exists()
        content = temp_file.read_text()
        assert content == "Column1, Column2, Column3\n"
        
    def test_append_rows_streaming(self, temp_file):
        """Test appending rows in streaming mode."""
        handler = OutputHandler(temp_file, mode='streaming')
        header = "Column1, Column2"
        handler.initialize_file(header)
        
        row1 = "7, 8"
        row2 = "9, 10"
        handler.append_row(row1)
        handler.append_row(row2)
        
        assert temp_file.exists()
        content = temp_file.read_text()
        assert content == "Column1, Column2\n7, 8\n9, 10\n"
        
    def test_append_rows_collect(self, temp_file):
        """Test appending rows in collect mode."""
        handler = OutputHandler(temp_file, mode='collect')
        header = "Column1, Column2"
        handler.initialize_file(header)
        
        row1 = "7, 8"
        row2 = "9, 10"
        handler.append_row(row1)
        handler.append_row(row2)
        
        assert temp_file.exists()
        content = temp_file.read_text()
        assert content == "Column1, Column2\n"
        
    def test_finalize_output_collect(self, temp_file):
        """Test finalizing output in collect mode."""
        handler = OutputHandler(temp_file, mode='collect')
        header = "Column1, Column2"
        
        row1 = "11, 12"
        row2 = "13, 14"
        handler.append_row(row1)
        handler.append_row(row2)
        
        handler.finalize(header)
        
        assert temp_file.exists()
        content = temp_file.read_text().strip().split('\n')
        assert content[0] == header
        assert content[1] == row1
        assert content[2] == row2
    
    def test_append_matrix_streaming(self, temp_file):
        """Test appending a matrix in streaming mode."""
        handler = OutputHandler(temp_file, mode='streaming')
        
        # Create a simple 3x3 matrix
        matrix = [
            [1, 2, 3],
            [4, 5, 6],
            [7, 8, 9]
        ]
        
        handler.append_matrix(matrix, frame_idx=0)
        
        # Check that frame-specific file was created
        matrix_file = temp_file.parent / "output_frame0.csv"
        assert matrix_file.exists()
        
        content = matrix_file.read_text()
        expected = "1,2,3\n4,5,6\n7,8,9\n"
        assert content == expected
    
    def test_append_matrix_collect(self, temp_file):
        """Test appending a matrix in collect mode."""
        handler = OutputHandler(temp_file, mode='collect')
        
        # Create two matrices for different frames
        matrix1 = [
            [1, 2],
            [3, 4]
        ]
        matrix2 = [
            [5, 6],
            [7, 8]
        ]
        
        handler.append_matrix(matrix1, frame_idx=0)
        handler.append_matrix(matrix2, frame_idx=1)
        
        # In collect mode, matrices should be stored, not written yet
        assert not temp_file.exists()
        assert len(handler.accumulated_matrices) == 2
        assert 0 in handler.accumulated_matrices
        assert 1 in handler.accumulated_matrices
        assert handler.accumulated_matrices[0] == matrix1
        assert handler.accumulated_matrices[1] == matrix2
    
    def test_finalize_matrices(self, temp_file):
        """Test finalizing matrices in collect mode."""
        handler = OutputHandler(temp_file, mode='collect')
        
        # Add multiple matrices
        matrix1 = [[1, 2], [3, 4]]
        matrix2 = [[5, 6], [7, 8]]
        matrix3 = [[9, 10], [11, 12]]
        
        handler.append_matrix(matrix1, frame_idx=0)
        handler.append_matrix(matrix2, frame_idx=5)
        handler.append_matrix(matrix3, frame_idx=10)
        
        # Finalize to write all matrices
        handler.finalize_matrices()
        
        # Check that frame-specific files were created
        frame0_file = temp_file.parent / "output_frame0.csv"
        frame5_file = temp_file.parent / "output_frame5.csv"
        frame10_file = temp_file.parent / "output_frame10.csv"
        
        assert frame0_file.exists()
        assert frame5_file.exists()
        assert frame10_file.exists()
        
        # Verify content
        assert frame0_file.read_text() == "1,2\n3,4\n"
        assert frame5_file.read_text() == "5,6\n7,8\n"
        assert frame10_file.read_text() == "9,10\n11,12\n"
    
    def test_append_matrix_streaming_multiple_frames(self, temp_file):
        """Test appending multiple matrices in streaming mode."""
        handler = OutputHandler(temp_file, mode='streaming')
        
        matrix1 = [[1, 2], [3, 4]]
        matrix2 = [[5, 6], [7, 8]]
        
        handler.append_matrix(matrix1, frame_idx=0)
        handler.append_matrix(matrix2, frame_idx=1)
        
        # Check that both frame-specific files were created immediately
        frame0_file = temp_file.parent / "output_frame0.csv"
        frame1_file = temp_file.parent / "output_frame1.csv"
        
        assert frame0_file.exists()
        assert frame1_file.exists()
        assert frame0_file.read_text() == "1,2\n3,4\n"
        assert frame1_file.read_text() == "5,6\n7,8\n"
 
# StructureAnalyzer tests # ------------------------------------------------------------       
class DummyAnalyzer(StructureAnalyzer):
    """A dummy analyzer for testing purposes."""
    
    def compute(self, graph):
        """Return simple test data."""
        return {"test": "data"}
    
    def render_output(self, data, frame_idx):
        """Write a simple row."""
        if self.output_handler:
            self.output_handler.append_row(f"{frame_idx},{data['test']}")

    def finalize_output(self, header):
        """Finalize output file."""
        super().finalize_output(header)
        
class TestStructureAnalyzer:
    """Test the StructureAnalyzer base class and DummyAnalyzer."""
    
    @pytest.fixture
    def temp_file(self):
        """Fixture that creates a temporary file path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir) / "output.csv"
            
    def test_dummy_analyzer_init(self, temp_file):
        """Test initialization of DummyAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='collect')
        analyzer = DummyAnalyzer(output_handler)
        
        assert analyzer.output_handler == output_handler
        assert analyzer.frame_number == 0
        
    def test_dummy_analyzer_compute(self, temp_file):
        """Test compute method of DummyAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='collect')
        analyzer = DummyAnalyzer(output_handler)
        
        result = analyzer.compute(None)
        assert result == {"test": "data"}
        
    def test_dummy_analyzer_render_output(self, temp_file):
        """Test render_output method of DummyAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='collect')
        analyzer = DummyAnalyzer(output_handler)
        header = "Column3, Column2"
        
        data = {"test": "data"}
        analyzer.render_output(data, frame_idx=5)
        analyzer.finalize_output(header)
        
        assert temp_file.exists()
        content = temp_file.read_text()
        assert content == "Column3, Column2\n5,data\n"
        
    def test_dummy_analyzer_finalize_output(self, temp_file):
        """Test finalize_output method of DummyAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='collect')
        analyzer = DummyAnalyzer(output_handler)
        header = "Frame, Data"
        
        data1 = {"test": "data1"}
        data2 = {"test": "data2"}
        analyzer.render_output(data1, frame_idx=1)
        analyzer.render_output(data2, frame_idx=4)
        analyzer.finalize_output(header)
        
        assert temp_file.exists()
        content = temp_file.read_text().strip().split('\n')
        assert content[0] == header
        assert content[1] == "1,data1"
        assert content[2] == "4,data2"
        
    def test_structure_analyzer_run(self):
        """Test the run method of StructureAnalyzer."""
        output_handler = OutputHandler(None, mode='collect')
        analyzer = DummyAnalyzer(output_handler)
        
        results = analyzer.run(graph=None, frame_idx=2)
        assert results == {"test": "data"}
        
# HBondsAnalyzer tests # ------------------------------------------------------------
class TestHBondsAnalyzer:
    """Test the HBondsAnalyzer class."""
    
    # SetUp fixtures -----------------------------------------
    @pytest.fixture
    def temp_file(self):
        """Fixture that creates a temporary file path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir) / "output.csv"
            
    @pytest.fixture
    def sample_file1(self):
        return "test_structures/Hbonds/01.xyz"
    
    @pytest.fixture
    def sample_file2(self):
        return "test_structures/Hbonds/02.xyz"
    
    @pytest.fixture
    def sample_file3(self):
        return "test_structures/Hbonds/03.xyz"
    
    @pytest.fixture
    def sample_file4(self):
        return "test_structures/Hbonds/04.xyz"
    
    @pytest.fixture
    def sample_file5(self):
        return "test_structures/Hbonds/05.xyz"
    
    @pytest.fixture
    def sample_file6(self):
        return "test_structures/Hbonds/06.xyz"
        
    # ---------------------------------------------------------
            
    def test_HBondsAnalyzer_init(self, temp_file):
        """Test initialization of HBondsAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='collect')
        cutoffs = [("O", 2.5, 3.5, 30), ("N", 2.7, 3.7, 25)]
        analyzer = HBondsAnalyzer(cutoffs, output_handler)
        
        assert analyzer.output_handler == output_handler
        assert analyzer.cutoffs == cutoffs
        
    # ---------------------------------------------------------
        
    def test_HBondsAnalyzer_compute(self, sample_file1):
        """Test compute function of HBondsAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file1))
        testGraph = graphs.GraphManager(xyz)
        
        cutoffs = [('O', 2.15, 3.8, 30)]
        analyzer = HBondsAnalyzer(cutoffs)
        hbonds = analyzer.compute(testGraph)
        assert len(hbonds) == 1
        assert hbonds[0][0] == 'O'
        assert hbonds[0][1] == 2.15
        assert hbonds[0][4] == 1
        
        cutoffs = [('O', 1.78, 3.5, 30)]
        analyzer = HBondsAnalyzer(cutoffs)
        hbonds = analyzer.compute(testGraph)
        assert len(hbonds) == 1
        assert hbonds[0][0] == 'O'
        assert hbonds[0][1] == 1.78
        assert hbonds[0][4] == 0
        
        cutoffs = [('O', 2.5, 3.8, 30)]
        analyzer = HBondsAnalyzer(cutoffs)
        hbonds = analyzer.compute(testGraph)
        assert len(hbonds) == 1
        assert hbonds[0][0] == 'O'
        assert hbonds[0][1] == 2.5
        assert hbonds[0][4] == 1

    def test2_HBondsAnalyzer_compute(self, sample_file2):
        xyz = next(readXYZ.readXYZ(sample_file2))
        testGraph = graphs.GraphManager(xyz)
        
        cutoffs = [('O', 2.15, 3.8, 30)]
        analyzer = HBondsAnalyzer(cutoffs)
        hbonds = analyzer.compute(testGraph)
        assert len(hbonds) == 1
        assert hbonds[0][0] == 'O'
        assert hbonds[0][1] == 2.15
        assert hbonds[0][4] == 6
        
    def test3_HBondsAnalyzer_compute(self, sample_file3):
        xyz = next(readXYZ.readXYZ(sample_file3))
        testGraph = graphs.GraphManager(xyz)
        
        cutoffs = [('O', 2.15, 3.8, 30)]
        analyzer = HBondsAnalyzer(cutoffs)
        hbonds = analyzer.compute(testGraph)
        assert len(hbonds) == 1
        assert hbonds[0][0] == 'O'
        assert hbonds[0][1] == 2.15
        assert hbonds[0][4] == 8
        
        cutoffs = [('N',2.5, 3.8, 30), ('O',2.15, 3.9, 30)]
        analyzer = HBondsAnalyzer(cutoffs)
        hbonds = analyzer.compute(testGraph)
        assert len(hbonds) == 2
        assert hbonds[0][0] == 'N'
        assert hbonds[0][1] == 2.5
        assert hbonds[0][4] == 0
        assert hbonds[1][0] == 'O'
        assert hbonds[1][1] == 2.15
        assert hbonds[1][4] == 8
        
        cutoffs = [('O', 2.15, 3.0, 30)]
        analyzer = HBondsAnalyzer(cutoffs)
        hbonds = analyzer.compute(testGraph)
        assert len(hbonds) == 1
        assert hbonds[0][0] == 'O'
        assert hbonds[0][1] == 2.15
        assert hbonds[0][4] == 7
        
        cutoffs = [('O', 2.15, 3.0, 31)]
        analyzer = HBondsAnalyzer(cutoffs)
        hbonds = analyzer.compute(testGraph)
        assert len(hbonds) == 1
        assert hbonds[0][0] == 'O'
        assert hbonds[0][1] == 2.15
        assert hbonds[0][4] == 8
        
        cutoffs = [('O', 2.15, 3.0, 40)]
        analyzer = HBondsAnalyzer(cutoffs)
        hbonds = analyzer.compute(testGraph)
        assert len(hbonds) == 1
        assert hbonds[0][0] == 'O'
        assert hbonds[0][1] == 2.15
        assert hbonds[0][4] == 9
        
        cutoffs = [('O', 2.0, 2.944, 30)]
        analyzer = HBondsAnalyzer(cutoffs)
        hbonds = analyzer.compute(testGraph)
        assert len(hbonds) == 1
        assert hbonds[0][0] == 'O'
        assert hbonds[0][1] == 2.0
        assert hbonds[0][2] == 2.944
        assert hbonds[0][3] == 30
        assert hbonds[0][4] == 5
        
    def test4_HBondsAnalyzer_compute(self, sample_file4):
        xyz = next(readXYZ.readXYZ(sample_file4))
        testGraph = graphs.GraphManager(xyz)
        
        cutoffs = [('N',2.1, 3.8, 30), ('O',2.15, 3.8, 30)]
        analyzer = HBondsAnalyzer(cutoffs)
        hbonds = analyzer.compute(testGraph)
        assert len(hbonds) == 2
        assert hbonds[0][0] == 'N'
        assert hbonds[0][1] == 2.1
        assert hbonds[0][4] == 2
        assert hbonds[1][0] == 'O'
        assert hbonds[1][1] == 2.15
        assert hbonds[1][4] == 0
        
    def test4_HBondsAnalyzer_compute(self, sample_file5):
        xyz = next(readXYZ.readXYZ(sample_file5))
        testGraph = graphs.GraphManager(xyz)
        
        cutoffs = [('N',2.1, 3.8, 30), ('O',2.15, 3.8, 30)]
        analyzer = HBondsAnalyzer(cutoffs)
        hbonds = analyzer.compute(testGraph)
        assert len(hbonds) == 2
        assert hbonds[0][0] == 'N'
        assert hbonds[0][1] == 2.1
        assert hbonds[0][4] == 8
        assert hbonds[1][0] == 'O'
        assert hbonds[1][1] == 2.15
        assert hbonds[1][4] == 0
        
    def test5_HBondsAnalyzer_compute(self, sample_file6):
        xyz = next(readXYZ.readXYZ(sample_file6))
        testGraph = graphs.GraphManager(xyz)
        
        cutoffs = [('O', 2.5, 3.8, 30)]
        analyzer = HBondsAnalyzer(cutoffs)
        hbonds = analyzer.compute(testGraph)
        assert len(hbonds) == 1
        assert hbonds[0][0] == 'O'
        assert hbonds[0][1] == 2.5
        assert hbonds[0][4] == 1
        
    # ---------------------------------------------------------
        
    def test_HBondsAnalyzer_render_finalize_output(self, temp_file):
        """Test render_output and render_output methods of HBondsAnalyzer."""
        # Mock the get_hbonds method of the GraphManager
        # -> we dont call the real graph methods here, since we only 
        #    want to test render_output
        mock_graph = Mock()
        mock_graph.get_hbonds.return_value = [(221, 3), (336, 129), (1, 18)]
        cutoffs = [('O', 2.5, 3.8, 30)]
        
        output_handler = OutputHandler(temp_file, mode='collect')
        analyzer = HBondsAnalyzer(cutoffs, output_handler)
        
        results = analyzer.compute(mock_graph)
        analyzer.render_output(results, frame_idx=18)
        
        results = analyzer.compute(mock_graph)
        analyzer.render_output(results, frame_idx=21)
        analyzer.finalize_output()
        
        assert temp_file.exists()
        content = temp_file.read_text()
        assert content == "Frame, Element Type, H-Acceptor dist / Å, Donor-Acceptor dist / Å, Angle cutoff / °, Number of Hydrogen Bonds\n18,O,2.500,3.800,30.000,3\n21,O,2.500,3.800,30.000,3\n"
        
# AnisotropyAnalyzer tests # ------------------------------------------------------------

class TestAnisotropyAnalyzer:
    """Test the AnisotropyAnalyzer class."""
    
    # SetUp fixtures -----------------------------------------
    @pytest.fixture
    def temp_file(self):
        """Fixture that creates a temporary file path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir) / "output.csv"
            
    @pytest.fixture
    def sample_file1(self):
        return "test_structures/Anisotropy/01.xyz"
    
    @pytest.fixture
    def sample_file2(self):
        return "test_structures/Anisotropy/02.xyz"
    
    @pytest.fixture
    def sample_file3(self):
        return "test_structures/Anisotropy/03.xyz"
    
    # --------------------------------------------------------
    
    def test_AnisotropyAnalyzer_init(self, temp_file):
        """Test initialization of HBondsAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='streaming')
        analyzer = AnisotropyAnalyzer(output_handler)
        
        assert analyzer.output_handler == output_handler
        
    # --------------------------------------------------------
    
    def test_AnisotropyAnalyzer_compute(self, temp_file):
        """Test compute function of AnisotropyAnalyzer."""
        output_handler = OutputHandler(temp_file)
        analyzer = AnisotropyAnalyzer(output_handler)
        
        mock_graph = Mock()
        
        with patch(
        "src.MakroLyzer.structure_modules.Anisotropy.get_Gtensor_eigVal_eigVec",
        return_value=(np.array([3.2, 3.5, 1.4]), None),
        ) as mock_helper:
            kappa_sq = analyzer.compute(mock_graph)
        mock_helper.assert_called_once_with(mock_graph)
        assert kappa_sq == pytest.approx(0.058985, rel=1e-5)
        
        with patch(
        "src.MakroLyzer.structure_modules.Anisotropy.get_Gtensor_eigVal_eigVec",
        return_value=(np.array([2, 2, 2]), None),
        ) as mock_helper:
            kappa_sq = analyzer.compute(mock_graph)
        mock_helper.assert_called_once_with(mock_graph)
        assert kappa_sq == pytest.approx(0.00000, rel=1e-5)
        
        with patch(
        "src.MakroLyzer.structure_modules.Anisotropy.get_Gtensor_eigVal_eigVec",
        return_value=(np.array([10, 0.1, 0.1]), None),
        ) as mock_helper:
            kappa_sq = analyzer.compute(mock_graph)
        mock_helper.assert_called_once_with(mock_graph)
        assert kappa_sq == pytest.approx(0.942042, rel=1e-5)
        
    # Linear Chain
    def test2_AnisotropyAnalyzer_compute(self, sample_file1):
        """Test compute function of AnisotropyAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file1))
        testGraph = graphs.GraphManager(xyz)
        
        analyzer = AnisotropyAnalyzer()
        kappa_sq = analyzer.compute(testGraph)
        assert kappa_sq == pytest.approx(1.0, abs=1e-3)
        
    # Planar - high symmetry
    def test3_AnisotropyAnalyzer_compute(self, sample_file2):
        """Test compute function of AnisotropyAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file2))
        testGraph = graphs.GraphManager(xyz)
        
        analyzer = AnisotropyAnalyzer()
        kappa_sq = analyzer.compute(testGraph)
        assert kappa_sq == pytest.approx(0.25, abs=1e-3)
        
    # 3D structure - high symmetry
    def test4_AnisotropyAnalyzer_compute(self, sample_file3):
        """Test compute function of AnisotropyAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file3))
        testGraph = graphs.GraphManager(xyz)
        
        analyzer = AnisotropyAnalyzer()
        kappa_sq = analyzer.compute(testGraph)
        assert kappa_sq == pytest.approx(0.00, abs=1e-3)
        
    # --------------------------------------------------------
    
    def test_AnisotropyAnalyzer_render_finalize_output(self, temp_file):
        """Test render_output and render_output methods of AnisotropyAnalyzer."""
        output_handler = OutputHandler(temp_file)
        analyzer = AnisotropyAnalyzer(output_handler)
        
        mock_graph = Mock()
        
        with patch(
        "src.MakroLyzer.structure_modules.Anisotropy.get_Gtensor_eigVal_eigVec",
        return_value=(np.array([3.2, 3.5, 1.4]), None),
        ) as mock_helper:
            kappa_sq = analyzer.compute(mock_graph)
        mock_helper.assert_called_once_with(mock_graph)
        
        analyzer.render_output(kappa_sq, frame_idx=12)
        
        with patch(
        "src.MakroLyzer.structure_modules.Anisotropy.get_Gtensor_eigVal_eigVec",
        return_value=(np.array([2, 2, 2]), None),
        ) as mock_helper:
            kappa_sq = analyzer.compute(mock_graph)
        mock_helper.assert_called_once_with(mock_graph)
        
        analyzer.render_output(kappa_sq, frame_idx=113)
        analyzer.finalize_output()
        
        assert temp_file.exists()
        content = temp_file.read_text()
        assert content == "Frame, Anisotropy Factor (κ²)\n12,0.059\n113,0.000\n"
        
# AsphericityAnalyzer tests # ------------------------------------------------------------

class TestAsphericityAnalyzer:
    """Test the AsphericityAnalyzer class."""
    
    # SetUp fixtures -----------------------------------------
    @pytest.fixture
    def temp_file(self):
        """Fixture that creates a temporary file path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir) / "hbond_output.csv"
            
    @pytest.fixture
    def sample_file1(self):
        return "test_structures/Asphericity/01.xyz"
    
    # --------------------------------------------------------
    
    def test_AsphericityAnalyzer_init(self, temp_file):
        """Test initialization of HBondsAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='streaming')
        analyzer = AsphericityAnalyzer(output_handler)
        
        assert analyzer.output_handler == output_handler
        
    # --------------------------------------------------------
    
    def test_AsphericityAnalyzer_compute(self, temp_file):
        """Test compute function of AsphericityAnalyzer."""
        output_handler = OutputHandler(temp_file)
        analyzer = AsphericityAnalyzer(output_handler)
        
        mock_graph = Mock()
        
        with patch(
        "src.MakroLyzer.structure_modules.Asphericity.get_Gtensor_eigVal_eigVec",
        return_value=(np.array([0, 102, 0]), None),
        ) as mock_helper:
            b = analyzer.compute(mock_graph)
        mock_helper.assert_called_once_with(mock_graph)
        assert b == pytest.approx(102.0000, rel=1e-5)
        
        with patch(
        "src.MakroLyzer.structure_modules.Asphericity.get_Gtensor_eigVal_eigVec",
        return_value=(np.array([2, 4, 1]), None),
        ) as mock_helper:
            b = analyzer.compute(mock_graph)
        mock_helper.assert_called_once_with(mock_graph)
        assert b == pytest.approx(2.50000, rel=1e-5)
        
        with patch(
        "src.MakroLyzer.structure_modules.Asphericity.get_Gtensor_eigVal_eigVec",
        return_value=(np.array([3.2, 3.2, 3.2]), None),
        ) as mock_helper:
            b = analyzer.compute(mock_graph)
        mock_helper.assert_called_once_with(mock_graph)
        assert b == pytest.approx(0.00000, rel=1e-5)
        
    def test2_AsphericityAnalyzer_compute(self, sample_file1):
        """Test compute function of AsphericityAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file1))
        testGraph = graphs.GraphManager(xyz)
        
        analyzer = AsphericityAnalyzer()
        b = analyzer.compute(testGraph)
        assert b == pytest.approx(0.0, abs=1e-3)
        
    # --------------------------------------------------------
    
    def test_AsphericityAnalyzer_render_finalize_output(self, temp_file):
        """Test render_output and render_output methods of AsphericityAnalyzer."""
        output_handler = OutputHandler(temp_file)
        analyzer = AsphericityAnalyzer(output_handler)
        
        mock_graph = Mock()
        
        with patch(
        "src.MakroLyzer.structure_modules.Asphericity.get_Gtensor_eigVal_eigVec",
        return_value=(np.array([0, 0, 100]), None),
        ) as mock_helper:
            b = analyzer.compute(mock_graph)
        mock_helper.assert_called_once_with(mock_graph)
        
        analyzer.render_output(b, frame_idx=12)
        
        with patch(
        "src.MakroLyzer.structure_modules.Asphericity.get_Gtensor_eigVal_eigVec",
        return_value=(np.array([4, 2, 1]), None),
        ) as mock_helper:
            b = analyzer.compute(mock_graph)
        mock_helper.assert_called_once_with(mock_graph)
        
        analyzer.render_output(b, frame_idx=113)
        analyzer.finalize_output()
        
        assert temp_file.exists()
        content = temp_file.read_text()
        assert content == "Frame, Asphericity Parameter (b)\n12,100.000\n113,2.500\n"
        
# RadiusOfGyrationAnalyzer tests # ----------------------------------------------------------

class TestRadiusOfGyrationAnalyzer:
    """Test the RadiusOfGyrationAnalyzer class."""

    # SetUp fixtures -----------------------------------------
    @pytest.fixture
    def temp_file(self):
        """Fixture that creates a temporary file path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir) / "output.csv"
            
    @pytest.fixture
    def sample_file1(self):
        return "test_structures/Rg/01.xyz"
    
    # ----------------------------------------------------------
    
    def test_RadiusOfGyrationAnalyzer_init(self, temp_file):
        output_handler = OutputHandler(temp_file, mode='streaming')
        analyzer = RadiusOfGyrationAnalyzer(output_handler)
        
        assert analyzer.output_handler == output_handler
        
    # ----------------------------------------------------------
    
    def test_RadiusOfGyrationAnalyzer_compute(self, temp_file):
        """Test compute function of RadiusOfGyrationAnalyzer."""
        output_handler = OutputHandler(temp_file)
        analyzer = RadiusOfGyrationAnalyzer(output_handler)
        
        mock_graph = Mock()
        mock_graph.get_all_coordinates.return_value = (
            None,
            np.array([(1.1, 22.0, 3.2), (1.4, 2.3, 19.0), (1.0, 2.0, 3.0)]),
        )
        mock_graph.get_com.return_value = (2.0, 4.2, 5.0)
        
        Rg = analyzer.compute(mock_graph)
        assert Rg == pytest.approx(13.300376, rel=1e-5) 
        
    def test2_RadiusOfGyrationAnalyzer_compute(self, sample_file1):
        """Test compute function of RadiusOfGyrationAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file1))
        testGraph = graphs.GraphManager(xyz)
        
        analyzer = RadiusOfGyrationAnalyzer()
        Rg = analyzer.compute(testGraph)
        assert Rg == pytest.approx(3.619, abs=1e-3)
    
    # ----------------------------------------------------------
    
    def test_RadiusOfGyrationAnalyzer_render_finalize_output(self, temp_file):
        """Test render_output and render_output methods of RadiusOfGyrationAnalyzer."""
        output_handler = OutputHandler(temp_file)
        analyzer = RadiusOfGyrationAnalyzer(output_handler)
        
        mock_graph = Mock()
        mock_graph.get_all_coordinates.return_value = (
            None,
            np.array([(1.1, 22.0, 3.2), (1.4, 2.3, 19.0), (1.0, 2.0, 3.0)]),
        )
        mock_graph.get_com.return_value = (2.0, 4.2, 5.0)
        
        Rg = analyzer.compute(mock_graph)
        analyzer.render_output(Rg, frame_idx=11)
        analyzer.finalize_output()
        
        assert temp_file.exists()
        content = temp_file.read_text()
        assert content == "Frame, Rg / Å\n11,13.300\n"
        
# test MoleculeCountAnalyzer tests # ----------------------------------------------------------

class TestMoleculeCountAnalyzer:
    """Test the MoleculeCountAnalyzer class."""

    # SetUp fixtures -----------------------------------------
    @pytest.fixture
    def temp_file(self):
        """Fixture that creates a temporary file path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir) / "molecule_count_output.csv"
            
    @pytest.fixture
    def sample_file1(self):
        return "test_structures/MoleculeCount/01.xyz"
    
    # ---------------------------------------------------------
    
    def test_MoleculeCountAnalyzer_init(self, temp_file):
        """Test initialization of MoleculeCountAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='streaming')
        analyzer = MoleculeCountAnalyzer(output_handler)
        
        assert analyzer.output_handler == output_handler
        
    # --------------------------------------------------------    
    
    def test_MoleculeCountAnalyzer_compute(self, sample_file1):
        """Test compute function of MoleculeCountAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file1))
        testGraph = graphs.GraphManager(xyz)
        
        analyzer = MoleculeCountAnalyzer()
        result = analyzer.compute(testGraph)
        assert result[0][0] == 12 # Total molecules
        assert result[0][1] == 5  # Cains
        assert result[0][2] == 7  # Rings
        
    def test2_MoleculeCountAnalyzer_compute(self, temp_file):
        """Test compute function of MoleculeCountAnalyzer."""
        # Create mock subgraphs
        mock_chain1 = Mock()
        mock_chain1.find_longest_path.return_value = [0, 1, 2, 3]
        mock_chain1.remove_1order.return_value = mock_chain1
        mock_chain1.subgraph.return_value.copy.return_value = mock_chain1
        mock_chain1.degree.side_effect = lambda n: 1 if n in [0, 3] else 2
        mock_chain1.nodes.return_value = [0, 1, 2, 3]
        
        mock_ring1 = Mock()
        mock_ring1.find_longest_path.return_value = [4, 5, 6, 7]
        mock_ring1.remove_1order.return_value = mock_ring1
        mock_ring1.subgraph.return_value.copy.return_value = mock_ring1
        mock_ring1.degree.side_effect = lambda n: 2  # All nodes have degree 2 (ring)
        mock_ring1.nodes.return_value = [4, 5, 6, 7]
        
        mock_chain2 = Mock()
        mock_chain2.find_longest_path.return_value = [8, 9, 10]
        mock_chain2.remove_1order.return_value = mock_chain2
        mock_chain2.subgraph.return_value.copy.return_value = mock_chain2
        mock_chain2.degree.side_effect = lambda n: 1 if n in [8, 10] else 2
        mock_chain2.nodes.return_value = [8, 9, 10]
        
        # Create mock graph
        mock_graph = Mock()
        mock_graph.get_subgraphs.return_value = [mock_chain1, mock_ring1, mock_chain2]
        
        # Create analyzer and compute
        output_handler = OutputHandler(temp_file, mode='streaming')
        analyzer = MoleculeCountAnalyzer(output_handler)
        result = analyzer.compute(mock_graph)
        
        # Assertions
        assert result[0][0] == 3  # Total molecules
        assert result[0][1] == 2  # Chains
        assert result[0][2] == 1  # Rings
        
        # Verify that methods were called
        mock_graph.get_subgraphs.assert_called()
        assert mock_chain1.find_longest_path.call_count == 1
        assert mock_ring1.find_longest_path.call_count == 1
        assert mock_chain2.find_longest_path.call_count == 1
        
    # --------------------------------------------------------   
    
    def test_MoleculeCountAnalyzer_render_finalize_output(self, temp_file):
        """Test render_output and render_output methods of MoleculeCountAnalyzer."""
        # Create mock subgraphs
        mock_chain1 = Mock()
        mock_chain1.find_longest_path.return_value = [0, 1, 2, 3]
        mock_chain1.remove_1order.return_value = mock_chain1
        mock_chain1.subgraph.return_value.copy.return_value = mock_chain1
        mock_chain1.degree.side_effect = lambda n: 1 if n in [0, 2] else 2
        mock_chain1.nodes.return_value = [0, 1, 2, 3]
        
        mock_ring1 = Mock()
        mock_ring1.find_longest_path.return_value = [4, 5, 6, 7]
        mock_ring1.remove_1order.return_value = mock_ring1
        mock_ring1.subgraph.return_value.copy.return_value = mock_ring1
        mock_ring1.degree.side_effect = lambda n: 2  # All nodes have degree 2 (ring)
        mock_ring1.nodes.return_value = [4, 5, 6, 7]
        
        
        # Create mock graph
        mock_graph = Mock()
        mock_graph.get_subgraphs.return_value = [mock_chain1, mock_ring1]
        
        # Create analyzer and compute
        output_handler = OutputHandler(temp_file, mode='streaming')
        analyzer = MoleculeCountAnalyzer(output_handler)
        analyzer.initialize_output()
        result = analyzer.compute(mock_graph)
        analyzer.render_output(result, frame_idx=190)
        analyzer.finalize_output()
        
        assert temp_file.exists()
        content = temp_file.read_text()
        assert content == "Frame, Molecule count, Chain count, Ring count\n190,2,1,1\n"
        
# EndToEndDistanceAnalyzer tests # ------------------------------------------------------------

class TestEndToEndDistanceAnalyzer:
    """Test the EndToEndDistanceAnalyzer class."""
    
    # SetUp fixtures -----------------------------------------
    @pytest.fixture
    def temp_file(self):
        """Fixture that creates a temporary file path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir) / "output.csv"
            
    @pytest.fixture
    def sample_file1(self):
        return "test_structures/EndToEnd/01.xyz"
    
    # --------------------------------------------------------
    
    def test_EndToEndDistanceAnalyzer_init(self, temp_file):
        """Test initializatio of EndToEndDistanceAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='streaming')
        analyzer = EndToEndDistanceAnalyzer(output_handler)
        
        assert analyzer.output_handler == output_handler
        
    # --------------------------------------------------------
        
    def test_EndToEndDistanceAnalyzer_compute(self, temp_file):
        """Test compute function of EndToEndDistanceAnalyzer."""
        mock_graph = Mock()
        mock_subgraph1 = Mock()
        mock_subgraph2 = Mock()
       
        mock_subgraph1.remove_1order.return_value = mock_subgraph1
        mock_subgraph2.remove_1order.return_value = mock_subgraph2
        mock_subgraph1.find_longest_path.return_value = [1,2,3]
        mock_subgraph2.find_longest_path.return_value = [1,2,3]
        mock_subgraph1.distance.return_value = 3.0
        mock_subgraph2.distance.return_value = 4.12
        
        mock_graph.get_subgraphs.return_value = [mock_subgraph1, mock_subgraph2]
        
        output_handler = OutputHandler(temp_file, mode='streaming')
        analyzer = EndToEndDistanceAnalyzer(output_handler)
        result = analyzer.compute(mock_graph)
        
        assert result[0] == 3.0
        assert result[1] == 4.12
        
    def test2_EndToEndDistanceAnalyzer_compute(self, temp_file):
        """Test compute function of EndToEndDistanceAnalyzer."""
        mock_graph = Mock()
        mock_subgraph1 = Mock()
        mock_subgraph2 = Mock()
       
        mock_subgraph1.remove_1order.return_value = mock_subgraph1
        mock_subgraph2.remove_1order.return_value = mock_subgraph2
        mock_subgraph1.find_longest_path.return_value = [1]
        mock_subgraph2.find_longest_path.return_value = [1,2,3]
        mock_subgraph1.distance.return_value = 3.0
        mock_subgraph2.distance.return_value = 4.12
        
        mock_graph.get_subgraphs.return_value = [mock_subgraph1, mock_subgraph2]
        
        output_handler = OutputHandler(temp_file, mode='streaming')
        analyzer = EndToEndDistanceAnalyzer(output_handler)
        result = analyzer.compute(mock_graph)
        
        assert result[0] == 4.12
        
    def test3_EndToEndDistanceAnalyzer_compute(self, sample_file1):
        """Test compute function of EndToEndDistanceAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file1))
        testGraph = graphs.GraphManager(xyz)
        
        analyzer = EndToEndDistanceAnalyzer()
        result = analyzer.compute(testGraph)
        assert result[0] == pytest.approx(5.595, abs=1e-3)
        
    # --------------------------------------------------------
    
    def test_EndToEndDistanceAnalyzer_render_finalize_output(self, temp_file):
        """Test render_output and render_output methods of EndToEndDistanceAnalyzer."""
        mock_graph = Mock()
        mock_subgraph1 = Mock()
        mock_subgraph2 = Mock()
       
        mock_subgraph1.remove_1order.return_value = mock_subgraph1
        mock_subgraph2.remove_1order.return_value = mock_subgraph2
        mock_subgraph1.find_longest_path.return_value = [1,2,3]
        mock_subgraph2.find_longest_path.return_value = [1,2,3]
        mock_subgraph1.distance.return_value = 3.0
        mock_subgraph2.distance.return_value = 4.6
        
        mock_graph.get_subgraphs.return_value = [mock_subgraph1, mock_subgraph2]
        
        output_handler = OutputHandler(temp_file, mode='streaming')
        analyzer = EndToEndDistanceAnalyzer(output_handler)
        analyzer.initialize_output()
        result = analyzer.compute(mock_graph)
        analyzer.render_output(result, frame_idx=12)
        analyzer.finalize_output()
        
        assert temp_file.exists()
        content = temp_file.read_text()
        assert content == "Frame, End-to-End Distance (subgraphs) / Å\n12,3.000,4.600\n"
        
# OrderParameterAnylyzer tests # ------------------------------------------------------------------

class TestOrderParameterAnalyzer:
    """Test the OrderParameterAnalyzer class."""
    
    # SetUp fixtures -----------------------------------------
    @pytest.fixture
    def temp_file(self):
        """Fixture that creates a temporary file path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir) / "output.csv"
            
    @pytest.fixture
    def sample_file1(self):
        return "test_structures/OrderParameter/01.xyz"
    
    @pytest.fixture
    def sample_file2(self):
        return "test_structures/OrderParameter/02.xyz"
    
    # ----------------------------------------------------------
    
    def test_OrderParameterAnalyzer_init(self, temp_file):
        output_handler = OutputHandler(temp_file, mode='streaming')
        BoxSize = (12.0, 13.4, 2.0)
        NoCellsPerDim = (1, 3, 4)
        MolecularVectorLength = 4
        analyzer = OrderParameterAnalyzer(BoxSize, NoCellsPerDim, MolecularVectorLength, output_handler)
        
        assert analyzer.BoxSize == BoxSize
        assert analyzer.NoCellsPerDim == NoCellsPerDim
        assert analyzer.MolecularVectorLength == MolecularVectorLength
        assert analyzer.output_handler == output_handler
        
    def test2_OrderParameterAnalyzer_init(self, temp_file):
        output_handler = OutputHandler(temp_file, mode='streaming')
        BoxSize = [12.0, 13.4]
        NoCellsPerDim = 1.3
        MolecularVectorLength = [4.2, 19]
        
        # we expect the WRONG_INPUT_TYPE_OP_ERROR since the inputs do not have the correct types
        expected_failure = pytest.raises(
            ValueError,
            match="Error: Wrong input type provided for the Order Parameter calculation. "
                  "Please ensure the input is of the correct type.",
        )
        with expected_failure:
            analyzer = OrderParameterAnalyzer(BoxSize, NoCellsPerDim, MolecularVectorLength, output_handler)
        
    # ----------------------------------------------------------
    
    def test_OrderParameterAnalyzer_get_backbone_vectors(self, temp_file):
        """Test get_backbone_vectors method of OrderParameterAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='streaming')
        BoxSize = (10.0, 10.0, 10.0)
        NoCellsPerDim = (1, 1, 1)
        MolecularVectorLength = 3
        analyzer = OrderParameterAnalyzer(BoxSize, NoCellsPerDim, MolecularVectorLength, output_handler)
        
        # Mock graph and its methods
        mock_graph = Mock()
        mock_new_graph = Mock()
        mock_subgraph = Mock()
        
        # Mock the chain of method calls
        mock_graph.remove_1order.return_value = mock_new_graph
        mock_new_graph.get_subgraphs.return_value = [mock_subgraph]
        
        # Mock longest path
        mock_subgraph.find_longest_path.return_value = [0, 1, 2, 3, 4]
        
        # Mock the get_vectors_and_positions_along_path method
        # Returns dict with midpoint (tuple) as key and vectors (list) as value
        mock_path_dict = {
            (1.0, 2.0, 3.0): [np.array([1.0, 0.0, 0.0])],
            (2.0, 3.0, 4.0): [np.array([0.0, 1.0, 0.0])]
        }
        mock_subgraph.get_vectors_and_positions_along_path.return_value = mock_path_dict
        
        # Call the method
        result = analyzer.get_backbone_vectors(mock_graph)
        
        # Verify the method calls
        mock_graph.remove_1order.assert_called_once()
        mock_new_graph.update_degree.assert_called_once()
        mock_new_graph.get_subgraphs.assert_called_once()
        mock_subgraph.find_longest_path.assert_called_once()
        mock_subgraph.get_vectors_and_positions_along_path.assert_called_once_with([0, 1, 2, 3, 4], MolecularVectorLength)
        
        # Verify the result structure
        assert isinstance(result, dict)
        assert (1.0, 2.0, 3.0) in result
        assert (2.0, 3.0, 4.0) in result
        assert len(result[(1.0, 2.0, 3.0)]) == 1
        assert np.array_equal(result[(1.0, 2.0, 3.0)][0], np.array([1.0, 0.0, 0.0]))
    
    # ----------------------------------------------------------
    
    def test_OrderParameterAnalyzer_get_backbone_vectors_in_cell(self, temp_file):
        """Test get_backbone_vectors_in_cell method of OrderParameterAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='streaming')
        BoxSize = (10.0, 10.0, 10.0)
        NoCellsPerDim = (1, 1, 1)
        MolecularVectorLength = 3
        analyzer = OrderParameterAnalyzer(BoxSize, NoCellsPerDim, MolecularVectorLength, output_handler)
        
        # Define a cell
        cell = [(0.0, 5.0), (0.0, 5.0), (0.0, 5.0)]
        
        # Create VecAndPos dictionary with midpoints and vectors
        VecAndPos = {
            (1.0, 2.0, 3.0): [np.array([1.0, 0.0, 0.0]), np.array([0.5, 0.5, 0.0])],
            (2.0, 3.0, 4.0): [np.array([0.0, 1.0, 0.0])],
            (6.0, 7.0, 8.0): [np.array([0.0, 0.0, 1.0])]  # This is outside the cell
        }
        
        # Create array of midpoints for the tree
        midpoints = np.array([
            [1.0, 2.0, 3.0],
            [2.0, 3.0, 4.0],
            [6.0, 7.0, 8.0]
        ])
        
        # Mock the cKDTree
        mock_tree = Mock()
        # The tree.query_ball_point should return indices of points that could be in the cell
        # Indices 0 and 1 are within the radius, index 2 is outside
        mock_tree.query_ball_point.return_value = [[0, 1, 2]]
        
        # Call the method
        result = analyzer.get_backbone_vectors_in_cell(cell, VecAndPos, midpoints, mock_tree)
        
        # Verify the tree was queried correctly
        mock_tree.query_ball_point.assert_called_once()
        
        # Verify the result: only vectors from midpoints inside the cell should be included
        # (1.0, 2.0, 3.0) and (2.0, 3.0, 4.0) are inside cell bounds
        # (6.0, 7.0, 8.0) is outside cell bounds
        assert isinstance(result, list)
        assert len(result) == 3  # 2 vectors from first midpoint + 1 from second
        assert np.array_equal(result[0], np.array([1.0, 0.0, 0.0]))
        assert np.array_equal(result[1], np.array([0.5, 0.5, 0.0]))
        assert np.array_equal(result[2], np.array([0.0, 1.0, 0.0]))
    
    # ----------------------------------------------------------
    
    def test_OrderParameterAnalyzer_Q_prime_tensor(self, temp_file):
        """Test Q_prime_tensor method of OrderParameterAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='streaming')
        BoxSize = (10.0, 10.0, 10.0)
        NoCellsPerDim = (1, 1, 1)
        MolecularVectorLength = 3
        analyzer = OrderParameterAnalyzer(BoxSize, NoCellsPerDim, MolecularVectorLength, output_handler)
        
        # Test with aligned vectors along x-axis
        vectors = [
            np.array([1.0, 0.0, 0.0]),
            np.array([2.0, 0.0, 0.0]),
            np.array([1.5, 0.0, 0.0])
        ]
        
        result = analyzer.Q_prime_tensor(vectors)
        
        # Verify result is a 3x3 matrix
        assert result.shape == (3, 3)
        
        # For perfectly aligned vectors along x-axis, Q_prime should be symmetric
        assert np.allclose(result, result.T)
        assert result[0][0] == pytest.approx(1.0, abs=1e-3)
        assert result[1][1] == pytest.approx(-0.5, abs=1e-3)
        assert result[2][2] == pytest.approx(-0.5, abs=1e-3)
        
        # Test with random vectors
        random_vectors = [
            np.array([1.0, 2.0, 3.0]),
            np.array([4.0, 10.0, 6.0]),
            np.array([7.0, 8.0, 9.0])
        ]
        
        result_random = analyzer.Q_prime_tensor(random_vectors)
        
        # Verify result is a 3x3 matrix
        assert result_random.shape == (3, 3)
        
        # Should be symmetric, trace should be close to zero
        assert np.allclose(result_random, result_random.T)  
        assert np.isclose(np.trace(result_random), 0.0)
        
        assert result_random[0][0] == pytest.approx(-0.28536548, abs=1e-3)  
        assert result_random[1][1] == pytest.approx(0.13675296, abs=1e-3)
        assert result_random[2][2] == pytest.approx(0.14861252, abs=1e-3)
        assert result_random[0][1] == pytest.approx(0.34733742, abs=1e-3)
        assert result_random[0][2] == pytest.approx(0.34846136, abs=1e-3)
        assert result_random[1][2] == pytest.approx(0.59722115, abs=1e-3)
        
        # Test with empty vectors list - should return NaN matrix
        empty_vectors = []
        result_empty = analyzer.Q_prime_tensor(empty_vectors)
        
        assert result_empty.shape == (3, 3)
        assert np.all(np.isnan(result_empty))
        
        # Test with orthogonal vectors (isotropic case)
        orthogonal_vectors = [
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0])
        ]
        
        result_orthogonal = analyzer.Q_prime_tensor(orthogonal_vectors)
        
        # For isotropic distribution, Q_prime should be close to zero matrix
        assert result_orthogonal.shape == (3, 3)
        # The trace should be zero for Q_prime tensor
        assert np.isclose(np.trace(result_orthogonal), 0.0)
    
    # ----------------------------------------------------------
    
    def test_OrderParameterAnalyzer_compute(self, sample_file1):
        """Test compute function of OrderParameterAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file1))
        testGraph = graphs.GraphManager(xyz, 3000)
        
        BoxSize = (3000, 3000, 3000)
        NoCellsPerDim = (1, 1, 1)
        MolecularVectorLength = 2
        
        analyzer = OrderParameterAnalyzer(BoxSize, NoCellsPerDim, MolecularVectorLength)
        S_star = analyzer.compute(testGraph)
        assert S_star == pytest.approx(1, abs=5e-2)
        
    def test2_OrderParameterAnalyzer_compute(self, sample_file2):
        """Test compute function of OrderParameterAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file2))
        testGraph = graphs.GraphManager(xyz, 3000)
        
        BoxSize = (3000, 3000, 3000)
        NoCellsPerDim = (1, 1, 1)
        MolecularVectorLength = 2
        
        analyzer = OrderParameterAnalyzer(BoxSize, NoCellsPerDim, MolecularVectorLength)
        S_star = analyzer.compute(testGraph)
        assert S_star == pytest.approx(0, abs=5e-2)
        
    def test3_OrderParameterAnalyzer_compute(self, temp_file):
        """Test compute function of OrderParameterAnalyzer."""
        # Create mock graph and its methods
        mock_graph = Mock()
        mock_new_graph = Mock()
        mock_subgraph = Mock()
        
        # Mock the chain of method calls
        mock_graph.remove_1order.return_value = mock_new_graph
        mock_new_graph.get_subgraphs.return_value = [mock_subgraph]
        
        # Mock longest path
        mock_subgraph.find_longest_path.return_value = [0, 1, 2, 3, 4]
        
        # Mock the get_vectors_and_positions_along_path method
        # Returns dict with midpoint (tuple) as key and vectors (list) as value
        mock_path_dict = {
            (1.0, 2.0, 3.0): [np.array([1.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])],
            (2.0, 3.0, 4.0): [np.array([1.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])]
        }
        mock_subgraph.get_vectors_and_positions_along_path.return_value = mock_path_dict
        
        BoxSize = (10.0, 10.0, 9.0)
        NoCellsPerDim = (1, 2, 1)
        MolecularVectorLength = 3
        
        analyzer = OrderParameterAnalyzer(BoxSize, NoCellsPerDim, MolecularVectorLength)
        S_star = analyzer.compute(mock_graph)
        
        assert S_star == pytest.approx(1.0, abs=1e-3)
        
        # Mock the get_vectors_and_positions_along_path method
        # Returns dict with midpoint (tuple) as key and vectors (list) as value
        mock_path_dict = {
            (1.0, 2.0, 3.0): [np.array([1.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])],
            (7.0, 3.0, 4.0): [np.array([1.0, 0.0, 0.0]), np.array([1.0, 0.0, 0.0])]
        }
        mock_subgraph.get_vectors_and_positions_along_path.return_value = mock_path_dict
        
        BoxSize = (10.0, 10.0, 9.0)
        NoCellsPerDim = (2, 1, 1)
        MolecularVectorLength = 3
        
        analyzer = OrderParameterAnalyzer(BoxSize, NoCellsPerDim, MolecularVectorLength)
        S_star = analyzer.compute(mock_graph)
        
        assert np.isnan(S_star)
        
    # ----------------------------------------------------------
    
    def test_OrderParameterAnalyzer_render_finalize_output(self, temp_file):
        """Test render_output and finalize_output methods of OrderParameterAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='collect')
        BoxSize = (10.0, 10.0, 10.0)
        NoCellsPerDim = (1, 1, 1)
        MolecularVectorLength = 3
        analyzer = OrderParameterAnalyzer(BoxSize, NoCellsPerDim, MolecularVectorLength, output_handler)
        
        # Mock graph
        mock_graph = Mock()
        
        # Mock the compute method to return a specific S* value
        with patch.object(analyzer, 'compute', return_value=0.856):
            S_star = analyzer.compute(mock_graph)
            analyzer.render_output(S_star, frame_idx=1)
        
        with patch.object(analyzer, 'compute', return_value=0.742):
            S_star = analyzer.compute(mock_graph)
            analyzer.render_output(S_star, frame_idx=5)
        
        analyzer.finalize_output()
        
        assert temp_file.exists()
        content = temp_file.read_text()
        assert content == "Frame, Order Parameter S*\n1,0.856\n5,0.742\n"
        
# RamachandranAnalyzer tests # -----------------------------------------------------------------------------

class TestRamachandranAnalyzer:
    """Test the RamachandranAnalyzer class."""
    
    # SetUp fixtures -----------------------------------------
    @pytest.fixture
    def temp_file(self):
        """Fixture that creates a temporary file path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir) / "output.csv"
            
    @pytest.fixture
    def sample_file1(self):
        return "test_structures/Ramachandran/alpha.xyz"
    
    @pytest.fixture
    def sample_file2(self):
        return "test_structures/Ramachandran/beta.xyz"
    
    # ----------------------------------------------------------
    
    def test_RamachandranAnalyzer_init(self, temp_file):
        """Test initialization of RamachandranAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='streaming')
        analyzer = RamachandranAnalyzer(output_handler)
        
        assert analyzer.output_handler == output_handler
        
    # --------------------------------------------------------
    
    def test_RamachandranAnalyzer_compute(self, sample_file1):
        """Test compute function of RamachandranAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file1))
        testGraph = graphs.GraphManager(xyz)
        
        analyzer = RamachandranAnalyzer()
        matrix = analyzer.compute(testGraph)
        
        # Verify matrix is 360x360
        assert len(matrix) == 360
        assert len(matrix[0]) == 360
        
        # Check some expected values 
        assert matrix[89][181] == 1
        assert matrix[112][139] == 1
        assert matrix[113][166] == 1
        assert matrix[116][139] == 1
        assert matrix[116][145] == 1
        assert matrix[117][152] == 1
        
    def test2_RamachandranAnalyzer_compute(self, sample_file2):
        """Test compute function of RamachandranAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file2))
        testGraph = graphs.GraphManager(xyz)
        
        analyzer = RamachandranAnalyzer()
        matrix = analyzer.compute(testGraph)
        
        # Verify matrix is 360x360
        assert len(matrix) == 360
        assert len(matrix[0]) == 360
        
        # Check some expected values 
        assert matrix[28][328] == 1
        assert matrix[34][267] == 1
        assert matrix[59][290] == 1
        assert matrix[81][343] == 1
        assert matrix[82][277] == 1
        assert matrix[265][166] == 1
    
    def test3_RamachandranAnalyzer_compute(self, temp_file):
        """Test compute function of RamachandranAnalyzer."""
        # Mock graph and subgraph
        mock_graph = Mock()
        mock_subgraph = Mock()
        
        # Mock get_subgraphs to return one subgraph
        mock_graph.get_subgraphs.return_value = [mock_subgraph]
        
        # Mock AminoAcidBackbone to return a backbone with enough atoms
        # Backbone pattern: N-Calpha-CarbonylC-N-Calpha-CarbonylC-N-...
        # Loop: i=2, checks i < len-4, increments by 3
        # For 2 iterations: i=2 (2<8), i=5 (5<8), i=8 (8<8 stops)
        # Need len(backbone) >= 10
        mock_subgraph.AminoAcidBackbone.return_value = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
        
        # Mock dihedral to return specific angles
        # First iteration (i=2): phi for atoms (2,3,4,5), psi for atoms (3,4,5,6)
        # Second iteration (i=5): phi for atoms (5,6,7,8), psi for atoms (6,7,8,9)
        mock_graph.dihedral.side_effect = [-45, 135, -30, 120]  # Two phi/psi pairs
        
        analyzer = RamachandranAnalyzer()
        matrix = analyzer.compute(mock_graph)
        
        # Verify matrix is 360x360
        assert len(matrix) == 360
        assert len(matrix[0]) == 360
        
        # Check that specific positions were incremented
        # First iteration: phi = -45 + 180 = 135, psi = 135 + 180 = 315
        assert matrix[135][315] == 1
        # Second iteration: phi = -30 + 180 = 150, psi = 120 + 180 = 300
        assert matrix[150][300] == 1
        
        # Verify the dihedral method was called correctly (2 iterations * 2 calls each)
        assert mock_graph.dihedral.call_count == 4
        
    def test_RamachandranAnalyzer_compute_no_backbone(self):
        """Test compute with subgraph that has insufficient backbone."""
        mock_graph = Mock()
        mock_subgraph = Mock()
        
        mock_graph.get_subgraphs.return_value = [mock_subgraph]
        # Backbone too short - only 5 atoms (need at least 7)
        mock_subgraph.AminoAcidBackbone.return_value = [0, 1, 2, 3, 4]
        
        analyzer = RamachandranAnalyzer()
        matrix = analyzer.compute(mock_graph)
        
        # Matrix should be all zeros since backbone is too short
        assert all(all(val == 0 for val in row) for row in matrix)
        
        # dihedral should not have been called
        mock_graph.dihedral.assert_not_called()
    
    def test_RamachandranAnalyzer__compute_multiple_subgraphs(self):
        """Test compute with multiple subgraphs."""
        mock_graph = Mock()
        mock_subgraph1 = Mock()
        mock_subgraph2 = Mock()
        
        mock_graph.get_subgraphs.return_value = [mock_subgraph1, mock_subgraph2]
        
        # Both subgraphs have valid backbones
        mock_subgraph1.AminoAcidBackbone.return_value = [0, 1, 2, 3, 4, 5, 6]
        mock_subgraph2.AminoAcidBackbone.return_value = [10, 11, 12, 13, 14, 15, 16]
        
        # Mock dihedrals for both subgraphs
        mock_graph.dihedral.side_effect = [45, 90, -60, -45]
        
        analyzer = RamachandranAnalyzer()
        matrix = analyzer.compute(mock_graph)
        
        # Verify both phi/psi pairs were added
        # First: phi=45+180=225, psi=90+180=270
        assert matrix[225][270] == 1
        # Second: phi=-60+180=120, psi=-45+180=135
        assert matrix[120][135] == 1
        
        assert mock_graph.dihedral.call_count == 4
        
    # --------------------------------------------------------
    
    def test_RamachandranAnalyzer_render_output(self, temp_file):
        """Test render_output method of RamachandranAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='collect')
        analyzer = RamachandranAnalyzer(output_handler)
        
        # Create a simple 3x3 matrix
        test_matrix = [
            [0, 1, 2],
            [3, 4, 5],
            [6, 7, 8]
        ]
        
        # Call render_output
        analyzer.render_output(test_matrix, frame_idx=0)
        
        # Verify matrix was accumulated
        assert 0 in output_handler.accumulated_matrices
        assert output_handler.accumulated_matrices[0] == test_matrix
        
        # No file should be written yet in collect mode
        assert not temp_file.exists()
    
    def test_RamachandranAnalyzer_finalize_output(self, temp_file):
        """Test finalize_output method of RamachandranAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='collect')
        analyzer = RamachandranAnalyzer(output_handler)
        
        # Add multiple matrices
        matrix1 = [[1, 2], [3, 4]]
        matrix2 = [[5, 6], [7, 8]]
        
        analyzer.render_output(matrix1, frame_idx=0)
        analyzer.render_output(matrix2, frame_idx=1)
        
        # Finalize should write all matrices
        analyzer.finalize_output()
        
        # Check that frame files were created
        frame0_file = temp_file.parent / "output_frame0.csv"
        frame1_file = temp_file.parent / "output_frame1.csv"
        
        assert frame0_file.exists()
        assert frame1_file.exists()
        
        # Verify content
        assert frame0_file.read_text() == "1,2\n3,4\n"
        assert frame1_file.read_text() == "5,6\n7,8\n"
        
# DihedralsAnalyzer tests # ------------------------------------------------------------------------
class TestDihedralsAnalyzer:
    """Test the DihedralsAnalyzer class."""
    
    # SetUp fixtures -----------------------------------------
    @pytest.fixture
    def temp_file(self):
        """Fixture that creates a temporary file path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir) / "output.csv"
            
    @pytest.fixture
    def sample_file1(self):
        return "test_structures/Dihedrals/01.xyz"
    
    @pytest.fixture
    def sample_file2(self):
        return "test_structures/Dihedrals/02.xyz"
    
    @pytest.fixture
    def sample_file3(self):
        return "test_structures/Dihedrals/03.xyz"
    
    @pytest.fixture
    def sample_file4(self):
        return "test_structures/Dihedrals/cis.xyz"
    
    @pytest.fixture
    def sample_file5(self):
        return "test_structures/Dihedrals/cis2.xyz"
    
    @pytest.fixture
    def sample_file6(self):
        return "test_structures/Dihedrals/trans.xyz"
    
    # ---------------------------------------------------------
    
    def test_DihedralsAnalyzer_init(self, temp_file):
        """Test initialization of DihedralsAnalyzer."""
        dihedral_output_handler = OutputHandler(temp_file, mode='collect')
        dihedral_list_output_handler = OutputHandler(temp_file, mode='streaming')
        cisTrans_output_handler = OutputHandler(temp_file, mode='streaming')
        sign = "nonabs"
        analyzer = DihedralsAnalyzer(
            dihedral_output_handler=dihedral_output_handler,
            dihedral_list_output_handler=dihedral_list_output_handler,
            cistrans_output_handler=cisTrans_output_handler,
            dihedral_range=sign
        )
        
        assert analyzer.dihedral_handler == dihedral_output_handler
        assert analyzer.dihedral_list_handler == dihedral_list_output_handler
        assert analyzer.cistrans_handler == cisTrans_output_handler
        assert analyzer.sign == True
        
    # ---------------------------------------------------------
    
    def test_DihedralsAnalyzer_get_all_dihedrals_mock(self):
        """Test _get_all_dihedrals method of DihedralsAnalyzer with mock data."""
        analyzer = DihedralsAnalyzer(dihedral_range='abs')
        
        # Create mock graph
        mock_graph = Mock()
        mock_prepared_graph = Mock()
        mock_subgraph1 = Mock()
        mock_subgraph2 = Mock()
        
        # Set up the graph preparation chain
        mock_graph.remove_1order.return_value = mock_prepared_graph
        mock_prepared_graph.get_subgraphs.return_value = [mock_subgraph1, mock_subgraph2]
        
        # Set up subgraph1: longest path with 4 nodes
        mock_subgraph1.find_longest_path.return_value = [0, 1, 2, 3]
        mock_subgraph1.dihedral.side_effect = [65.4, 75.0]  # Two dihedrals
        
        # Set up subgraph2: longest path with 5 nodes
        mock_subgraph2.find_longest_path.return_value = [0, 1, 2, 3, 4]
        mock_subgraph2.dihedral.side_effect = [44.8, 55.0]  # Two dihedrals
        
        dihedrals, dihedral_list = analyzer._get_all_dihedrals(mock_graph)
        
        # Verify structure
        assert isinstance(dihedrals, list)
        assert isinstance(dihedral_list, list)
        assert len(dihedral_list) == 2  # Two subgraphs
        
        # Verify angles per subgraph
        assert dihedral_list[0] == [65] # 75 is not included since 4 nodes -> 1 dihedral
        assert dihedral_list[1] == [45, 55]
        
        # Verify angle counts
        angle_dict = {angle: count for angle, count in dihedrals}
        assert angle_dict[45] == 1
        assert angle_dict[55] == 1
        assert angle_dict[65] == 1  
    
    def test_DihedralsAnalyzer_get_all_dihedrals(self, sample_file1):
        """Test _get_all_dihedrals method of DihedralsAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file1))
        testGraph = graphs.GraphManager(xyz)
        
        analyzer = DihedralsAnalyzer(dihedral_range='abs')
        dihedrals, dihedral_list = analyzer._get_all_dihedrals(testGraph)
        
        # Verify structure: list of (angle, count) tuples
        assert isinstance(dihedrals, list)
        assert all(isinstance(item, tuple) and len(item) == 2 for item in dihedrals)
        
        # Verify angles are in correct range (0-180)
        for angle, count in dihedrals:
            assert 0 <= angle <= 180
            assert count >= 0
        
        # Verify dihedral_list structure: list of lists per subgraph
        assert isinstance(dihedral_list, list)
        assert all(isinstance(subgraph_angles, list) for subgraph_angles in dihedral_list)
        
        # Verify specific angles match expected values
        angle_dict = {angle: count for angle, count in dihedrals}
        assert angle_dict[65] == 1
        assert angle_dict[177] == 2
        assert angle_dict[180] == 1
        
    def test_DihedralsAnalyzer_get_cistrans_from_dihedrals_mock(self):
        """Test _get_cistrans_from_dihedrals method of DihedralsAnalyzer with mock data."""
        analyzer = DihedralsAnalyzer(dihedral_range='abs')
        
        # Create mock dihedral list with known angles
        # 0-90 = cis, 90-180 = trans
        dihedral_list = [
            [30, 45, 90],      # 2 cis, 1 on boundary
            [100, 120, 150],   # 3 trans
            [60, 85]           # 2 cis
        ]
        
        cistrans = analyzer._get_cistrans_from_dihedrals(dihedral_list)
        
        # Verify structure
        assert isinstance(cistrans, list)
        assert len(cistrans) == 2
        assert cistrans[0][0] == 'Cis'
        assert cistrans[1][0] == 'Trans'
        
        # Verify counts: 4 cis (30, 45, 60, 85) and 3 trans (100, 120, 150), 90 is included in cis (0 <= d <= 90)
        assert cistrans[0][1] == 5  # Cis count
        assert cistrans[1][1] == 3  # Trans count
        
    def test_DihedralsAnalyzer_get_cistrans_from_dihedrals(self, sample_file1):
        """Test _get_cistrans_from_dihedrals method of DihedralsAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file1))
        testGraph = graphs.GraphManager(xyz)
        
        analyzer = DihedralsAnalyzer(dihedral_range='abs')
        _, dihedral_list = analyzer._get_all_dihedrals(testGraph)
        
        cistrans = analyzer._get_cistrans_from_dihedrals(dihedral_list)
        
        # Verify structure: list of (label, count) tuples
        assert isinstance(cistrans, list)
        assert len(cistrans) == 2
        assert cistrans[0][0] == 'Cis'
        assert cistrans[1][0] == 'Trans'
        
        # Verify counts match expected values
        # dihedral_list[0] = [180, 177, 65, 177] -> 1 cis (65), 3 trans (180, 177, 177)
        assert cistrans[0][1] == 1  # Cis count
        assert cistrans[1][1] == 3  # Trans count
        
    # ---------------------------------------------------------
    
    def test_DihedralsAnalyzer_compute_mock(self, temp_file):
        """Test compute function of DihedralsAnalyzer with mock data."""
        dihedral_handler = OutputHandler(temp_file, mode='collect')
        analyzer = DihedralsAnalyzer(
            dihedral_output_handler=dihedral_handler,
            dihedral_range='abs'
        )
        
        # Create mock graph
        mock_graph = Mock()
        mock_prepared_graph = Mock()
        mock_subgraph = Mock()
        
        # Set up the graph preparation chain
        mock_graph.remove_1order.return_value = mock_prepared_graph
        mock_prepared_graph.get_subgraphs.return_value = [mock_subgraph]
        
        # Set up subgraph: longest path with 5 nodes
        mock_subgraph.find_longest_path.return_value = [0, 1, 2, 3, 4]
        mock_subgraph.dihedral.side_effect = [50.0, 120.0]  # One cis, one trans
        
        result = analyzer.compute(mock_graph)
        
        # Verify result structure
        assert 'dihedrals' in result
        assert 'dihedral_list' in result
        assert 'cistrans' in result
        
        # Verify dihedrals
        dihedrals = result['dihedrals']
        assert isinstance(dihedrals, list)
        angle_dict = {angle: count for angle, count in dihedrals}
        assert angle_dict[50] == 1
        assert angle_dict[120] == 1
        
        # Verify dihedral_list
        dihedral_list = result['dihedral_list']
        assert dihedral_list[0] == [50, 120]
        
        # Verify cistrans: 50 is cis (0-90), 120 is trans (90-180)
        cistrans = result['cistrans']
        assert cistrans[0][0] == 'Cis'
        assert cistrans[0][1] == 1
        assert cistrans[1][0] == 'Trans'
        assert cistrans[1][1] == 1
        
    def test_DihedralsAnalyzer_compute(self, sample_file1):
        """Test compute function of DihedralsAnalyzer with real data."""
        xyz = next(readXYZ.readXYZ(sample_file1))
        testGraph = graphs.GraphManager(xyz)
        
        dihedral_range = "abs"
        analyzer = DihedralsAnalyzer(dihedral_range=dihedral_range)
        result = analyzer.compute(testGraph)
        
        # Unpack the dictionary
        dihedrals = result['dihedrals']
        dihedral_list = result['dihedral_list']
        cis_trans = result['cistrans']
        
        # Test dihedral counts
        for angle, count in dihedrals:
            if angle == 65:
                assert count == 1
            elif angle == 177:
                assert count == 2
            elif angle == 180:
                assert count == 1
            else:
                assert count == 0
                
        # test dihedral list
        assert len(dihedral_list[0]) == 4
        assert dihedral_list[0][0] == 180
        assert dihedral_list[0][1] == 177
        assert dihedral_list[0][2] == 65
        assert dihedral_list[0][3] == 177
                
        # test cis trans counts
        assert cis_trans[0][1] == 1
        assert cis_trans[1][1] == 3
        
        dihedral_range = "nonabs"
        analyzer = DihedralsAnalyzer(dihedral_range=dihedral_range)
        result = analyzer.compute(testGraph)
        
        # Unpack the dictionary
        dihedrals = result['dihedrals']
        dihedral_list = result['dihedral_list']
        cis_trans = result['cistrans']
        
        # Test dihedral counts
        for angle, count in dihedrals:
            if angle == -65:
                assert count == 1
            elif angle == -177:
                assert count == 2
            elif angle == 180:
                assert count == 1
            else: 
                assert count == 0
                
        # test dihedral list
        assert len(dihedral_list[0]) == 4
        assert dihedral_list[0][0] == 180
        assert dihedral_list[0][1] == -177
        assert dihedral_list[0][2] == -65
        assert dihedral_list[0][3] == -177
                
        # test cis trans counts
        assert cis_trans[0][1] == 1
        assert cis_trans[1][1] == 3
        
    def test2_DihedralsAnalyzer_compute(self, sample_file4):
        """Test compute function of DihedralsAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file4))
        testGraph = graphs.GraphManager(xyz)
        
        dihedral_range = "nonabs"
        analyzer = DihedralsAnalyzer(dihedral_range=dihedral_range)
        result = analyzer.compute(testGraph)
        
        # Unpack the dictionary
        dihedrals = result['dihedrals']
        dihedral_list = result['dihedral_list']
        cis_trans = result['cistrans']
        
        # Test dihedral counts
        for angle, count in dihedrals:
            if angle == -64:
                assert count == 1
            else:
                assert count == 0
                
        # test dihedral list
        assert len(dihedral_list[0]) == 1
        assert dihedral_list[0][0] == -64
                
        # test cis trans counts
        assert cis_trans[0][1] == 1
        assert cis_trans[1][1] == 0
        
    def test3_DihedralsAnalyzer_compute(self, sample_file5):
        """Test compute function of DihedralsAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file5))
        testGraph = graphs.GraphManager(xyz)
        
        dihedral_range = "nonabs"
        analyzer = DihedralsAnalyzer(dihedral_range=dihedral_range)
        result = analyzer.compute(testGraph)
        
        # Unpack the dictionary
        dihedrals = result['dihedrals']
        dihedral_list = result['dihedral_list']
        cis_trans = result['cistrans']
        
        # Test dihedral counts
        for angle, count in dihedrals:
            if angle == 4:
                assert count == 1
            else:
                assert count == 0
                
        # test dihedral list
        assert len(dihedral_list[0]) == 1
        assert dihedral_list[0][0] == 4
                
        # test cis trans counts
        assert cis_trans[0][1] == 1
        assert cis_trans[1][1] == 0
        
    def test4_DihedralsAnalyzer_compute(self, sample_file6):
        """Test compute function of DihedralsAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file6))
        testGraph = graphs.GraphManager(xyz)
        
        dihedral_range = "nonabs"
        analyzer = DihedralsAnalyzer(dihedral_range=dihedral_range)
        result = analyzer.compute(testGraph)
        
        # Unpack the dictionary
        dihedrals = result['dihedrals']
        dihedral_list = result['dihedral_list']
        cis_trans = result['cistrans']
        
        # Test dihedral counts
        for angle, count in dihedrals:
            if angle == 180:
                assert count == 1
            else:
                assert count == 0
                
        # test dihedral list
        assert len(dihedral_list[0]) == 1
        assert dihedral_list[0][0] == 180
                
        # test cis trans counts
        assert cis_trans[0][1] == 0
        assert cis_trans[1][1] == 1
        
    def test5_DihedralsAnalyzer_compute(self, sample_file2):
        """Test compute function of DihedralsAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file2))
        testGraph = graphs.GraphManager(xyz)
        
        dihedral_range = "abs"
        analyzer = DihedralsAnalyzer(dihedral_range=dihedral_range)
        result = analyzer.compute(testGraph)
        
        # Unpack the dictionary
        dihedrals = result['dihedrals']
        dihedral_list = result['dihedral_list']
        cis_trans = result['cistrans']
        
        # Test dihedral counts
        for angle, count in dihedrals:
            if angle == 2:
                assert count == 1
            elif angle == 64:
                assert count == 1
            elif angle == 66:
                assert count == 1
            elif angle == 81:
                assert count == 1
            elif angle == 172:
                assert count == 1
            elif angle == 175:
                assert count == 1
            elif angle == 177:
                assert count == 2
            elif angle == 179:
                assert count == 1
            elif angle == 180:
                assert count == 4
            else:
                assert count == 0
                
        # test cis trans counts
        assert cis_trans[0][1] == 4
        assert cis_trans[1][1] == 9
        
        dihedral_range = "nonabs"
        analyzer = DihedralsAnalyzer(dihedral_range=dihedral_range)
        result = analyzer.compute(testGraph)
        
        # Unpack the dictionary
        dihedrals = result['dihedrals']
        dihedral_list = result['dihedral_list']
        cis_trans = result['cistrans']
        
        # Test dihedral counts
        for angle, count in dihedrals:
            if angle == -175:
                assert count == 1
            elif angle == -81:
                assert count == 1
            elif angle == -2:
                assert count == 1
            elif angle == 172:
                assert count == 1
            elif angle == 64:
                assert count == 1
            elif angle == 66:
                assert count == 1
            elif angle == 177:
                assert count == 2
            elif angle == 179:
                assert count == 1
            elif angle == 180:
                assert count == 4
            else:
                assert count == 0
                
        # test cis trans counts
        assert cis_trans[0][1] == 4
        assert cis_trans[1][1] == 9
        
    def test6_DihedralsAnalyzer_compute(self, sample_file3):
        """Test compute function of DihedralsAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file3))
        testGraph = graphs.GraphManager(xyz)
        
        dihedral_range = "abs"
        analyzer = DihedralsAnalyzer(dihedral_range=dihedral_range)
        result = analyzer.compute(testGraph)
        
        # Unpack the dictionary
        dihedrals = result['dihedrals']
        dihedral_list = result['dihedral_list']
        cis_trans = result['cistrans']
        
        # Test dihedral counts
        for angle, count in dihedrals:
            if angle == 0:
                assert count == 1
            else:
                assert count == 0
                
        # test cis trans counts
        assert cis_trans[0][1] == 1
        assert cis_trans[1][1] == 0
        
    # ---------------------------------------------------------
        
    def test_DihedralsAnalyzer_render_finalize_output(self, temp_file):
        """Test render_output and finalize_output methods of DihedralsAnalyzer."""
        dihedral_handler = OutputHandler(temp_file, mode='collect')
        dihedral_list_handler = OutputHandler(temp_file.parent / "dihedral_list.csv", mode='collect')
        cistrans_handler = OutputHandler(temp_file.parent / "cistrans.csv", mode='collect')
        
        analyzer = DihedralsAnalyzer(
            dihedral_output_handler=dihedral_handler,
            dihedral_list_output_handler=dihedral_list_handler,
            cistrans_output_handler=cistrans_handler,
            dihedral_range='abs'
        )
        analyzer.initialize_output()
        
        # Create mock data and render it
        mock_graph = Mock()
        mock_prepared_graph = Mock()
        mock_subgraph = Mock()
        
        mock_graph.remove_1order.return_value = mock_prepared_graph
        mock_prepared_graph.get_subgraphs.return_value = [mock_subgraph]
        mock_subgraph.find_longest_path.return_value = [0, 1, 2, 3]
        mock_subgraph.dihedral.side_effect = [65.0, 75.0]
        
        # Render first frame
        result1 = analyzer.compute(mock_graph)
        analyzer.render_output(result1, frame_idx=0)
        
        # Reset mock for second frame
        mock_subgraph.dihedral.side_effect = [45.0, 85.0]
        result2 = analyzer.compute(mock_graph)
        analyzer.render_output(result2, frame_idx=1)
        
        # Finalize output
        analyzer.finalize_output()
        
        # Verify dihedral counts file
        assert temp_file.exists()
        dihedral_content = temp_file.read_text()
        lines = dihedral_content.strip().split('\n')
        assert 'Angle' in lines[0]
        assert 'Frame 0' in lines[0]
        assert 'Frame 1' in lines[0]
        
# ChemicalFormulaAnalyzer tests # ---------------------------------------------------------

class TestChemicalFormulaAnalyzer:
    """Test the ChemicalFormulaAnalyzer class."""
    
    # SetUp fixtures -----------------------------------------
    @pytest.fixture
    def temp_file(self):
        """Fixture that creates a temporary file path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir) / "output.csv"
            
    @pytest.fixture
    def sample_file1(self):
        return "test_structures/ChemicalFormula/01.xyz"
    
    # --------------------------------------------------------
    
    def test_ChemicalFormulaAnalyzer_init(self, temp_file):
        """Test initialization of ChemicalFormulaAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='streaming')
        analyzer = ChemicalFormulaAnalyzer(output_handler)
        
        assert analyzer.output_handler == output_handler
        
    # --------------------------------------------------------
    
    def test_ChemicalFormulaAnalyzer_compute(self, sample_file1):
        """Test compute function of ChemicalFormulaAnalyzer."""
        xyz = next(readXYZ.readXYZ(sample_file1))
        testGraph = graphs.GraphManager(xyz)
        formula = [('C1H6Mg1N1O1P1', 2), ('C1H5Mg1N1O1P1', 1)]
        
        analyzer = ChemicalFormulaAnalyzer()
        result = analyzer.compute(testGraph)
        for i in range(len(formula)):
            assert result[i][0] == formula[i][0]
            assert result[i][1] == formula[i][1]
    
    def test2_ChemicalFormulaAnalyzer_compute(self, temp_file):
        """Test compute function of ChemicalFormulaAnalyzer with mock graph."""
        # Create mock subgraphs with different chemical formulas
        mock_subgraph1 = Mock()
        mock_subgraph1.nodes = {0: {'element': 'C'}, 1: {'element': 'H'}, 2: {'element': 'O'}}
        
        mock_subgraph2 = Mock()
        mock_subgraph2.nodes = {3: {'element': 'C'}, 4: {'element': 'H'}}
        
        mock_subgraph3 = Mock()
        mock_subgraph3.nodes = {5: {'element': 'C'}, 6: {'element': 'H'}}
        
        mock_subgraph4 = Mock()
        mock_subgraph4.nodes = {5: {'element': 'C'}, 6: {'element': 'H'}, 7: {'element': 'H'}}
        
        # Create mock graph
        mock_graph = Mock()
        mock_graph.get_subgraphs.return_value = [mock_subgraph1, mock_subgraph2, mock_subgraph3, mock_subgraph4]
        
        # Create analyzer and compute
        output_handler = OutputHandler(temp_file, mode='streaming')
        analyzer = ChemicalFormulaAnalyzer(output_handler)
        result = analyzer.compute(mock_graph)
        
        # Assertions - result should be sorted by count (descending)
        # CH1O1 appears once, C1H1 appears twice
        assert result[0][0] == 'C1H1'  # Most frequent formula
        assert result[0][1] == 2
        assert result[1][0] == 'C1H1O1'  # Less frequent formula
        assert result[1][1] == 1
        assert result[2][0] == 'C1H2'
        assert result[2][1] == 1
        
        # Verify that get_subgraphs was called
        mock_graph.get_subgraphs.assert_called()
    
    def test_ChemicalFormulaAnalyzer_render_finalize_output(self, temp_file):
        """Test render_output and finalize_output methods of ChemicalFormulaAnalyzer."""
        # Create mock subgraphs with repeated formulas to test sorting
        mock_subgraph1 = Mock()
        mock_subgraph1.nodes = {0: {'element': 'C'}, 1: {'element': 'H'}}
        
        mock_subgraph2 = Mock()
        mock_subgraph2.nodes = {2: {'element': 'C'}, 3: {'element': 'H'}}
        
        mock_subgraph3 = Mock()
        mock_subgraph3.nodes = {4: {'element': 'C'}, 5: {'element': 'H'}, 6: {'element': 'O'}}
        
        # Create mock graph
        mock_graph = Mock()
        mock_graph.get_subgraphs.return_value = [mock_subgraph1, mock_subgraph2, mock_subgraph3]
        
        # Create analyzer and compute
        output_handler = OutputHandler(temp_file, mode='streaming')
        analyzer = ChemicalFormulaAnalyzer(output_handler)
        analyzer.initialize_output()
        result = analyzer.compute(mock_graph)
        analyzer.render_output(result, frame_idx=5)
        analyzer.finalize_output()
        
        assert temp_file.exists()
        content = temp_file.read_text()
        assert content == "Frame, Chemical Formula, Count\n5,C1H1,2\n5,C1H1O1,1\n"
    