import pytest
import tempfile
from pathlib import Path
import numpy as np

from unittest.mock import Mock, patch

from src.MakroLyzer.input_handling import readXYZ
from src.MakroLyzer.structure_modules import graphs

from src.MakroLyzer.structure_modules.structureBase import OutputHandler, StructureAnalyzer
from src.MakroLyzer.structure_modules.Hbonds import HBondsAnalyzer
from src.MakroLyzer.structure_modules.Anisotropy import AnisotropyAnalyzer
from src.MakroLyzer.structure_modules.Asphericity import AsphericityAnalyzer
from src.MakroLyzer.structure_modules.RadiusOfGyration import RadiusOfGyrationAnalyzer
from src.MakroLyzer.structure_modules.MoleculeCount import MoleculeCountAnalyzer

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
            yield Path(tmpdir) / "hbond_output.csv"
            
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
        """Test render_output and render_output methods of HBondsAnalyzer"""
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
            yield Path(tmpdir) / "hbond_output.csv"
            
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
        """Test compute function of AnisotropyAnalyzer"""
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
        """Test compute function of AnisotropyAnalyzer"""
        xyz = next(readXYZ.readXYZ(sample_file1))
        testGraph = graphs.GraphManager(xyz)
        
        analyzer = AnisotropyAnalyzer()
        kappa_sq = analyzer.compute(testGraph)
        assert kappa_sq == pytest.approx(1.0, abs=1e-3)
        
    # Planar - high symmetry
    def test3_AnisotropyAnalyzer_compute(self, sample_file2):
        """Test compute function of AnisotropyAnalyzer"""
        xyz = next(readXYZ.readXYZ(sample_file2))
        testGraph = graphs.GraphManager(xyz)
        
        analyzer = AnisotropyAnalyzer()
        kappa_sq = analyzer.compute(testGraph)
        assert kappa_sq == pytest.approx(0.25, abs=1e-3)
        
    # 3D structure - high symmetry
    def test4_AnisotropyAnalyzer_compute(self, sample_file3):
        """Test compute function of AnisotropyAnalyzer"""
        xyz = next(readXYZ.readXYZ(sample_file3))
        testGraph = graphs.GraphManager(xyz)
        
        analyzer = AnisotropyAnalyzer()
        kappa_sq = analyzer.compute(testGraph)
        assert kappa_sq == pytest.approx(0.00, abs=1e-3)
        
    # --------------------------------------------------------
    
    def test_AnisotropyAnalyzer_render_finalize_output(self, temp_file):
        """Test render_output and render_output methods of AnisotropyAnalyzer"""
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
        """Test compute function of AsphericityAnalyzer"""
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
        """Test compute function of AsphericityAnalyzer"""
        xyz = next(readXYZ.readXYZ(sample_file1))
        testGraph = graphs.GraphManager(xyz)
        
        analyzer = AsphericityAnalyzer()
        b = analyzer.compute(testGraph)
        assert b == pytest.approx(0.0, abs=1e-3)
        
    # --------------------------------------------------------
    
    def test_AsphericityAnalyzer_render_finalize_output(self, temp_file):
        """Test render_output and render_output methods of AsphericityAnalyzer"""
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
            yield Path(tmpdir) / "hbond_output.csv"
            
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
        """Test compute function of RadiusOfGyrationAnalyzer"""
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
        """Test compute function of RadiusOfGyrationAnalyzer"""
        xyz = next(readXYZ.readXYZ(sample_file1))
        testGraph = graphs.GraphManager(xyz)
        
        analyzer = RadiusOfGyrationAnalyzer()
        Rg = analyzer.compute(testGraph)
        assert Rg == pytest.approx(3.619, abs=1e-3)
    
    # ----------------------------------------------------------
    
    def test_RadiusOfGyrationAnalyzer_render_finalize_output(self, temp_file):
        """Test render_output and render_output methods of RadiusOfGyrationAnalyzer"""
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