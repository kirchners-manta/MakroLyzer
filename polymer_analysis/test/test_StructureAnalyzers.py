import pytest
import tempfile
from pathlib import Path

from unittest.mock import Mock, MagicMock

from src.MakroLyzer.input_handling import readXYZ
from src.MakroLyzer.structure_modules import graphs

from src.MakroLyzer.structure_modules.structureBase import OutputHandler, StructureAnalyzer
from src.MakroLyzer.structure_modules.Hbonds import HBondsAnalyzer



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
            
    def test_hbonds_analyzer_init(self, temp_file):
        """Test initialization of HBondsAnalyzer."""
        output_handler = OutputHandler(temp_file, mode='collect')
        cutoffs = [("O", 2.5, 3.5, 30), ("N", 2.7, 3.7, 25)]
        analyzer = HBondsAnalyzer(cutoffs, output_handler)
        
        assert analyzer.output_handler == output_handler
        assert analyzer.cutoffs == cutoffs
        
    # ---------------------------------------------------------
        
    def test_hbonds_analyzer_compute(self, sample_file1):
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

    def test2_hbonds_analyzer_compute(self, sample_file2):
        xyz = next(readXYZ.readXYZ(sample_file2))
        testGraph = graphs.GraphManager(xyz)
        
        cutoffs = [('O', 2.15, 3.8, 30)]
        analyzer = HBondsAnalyzer(cutoffs)
        hbonds = analyzer.compute(testGraph)
        assert len(hbonds) == 1
        assert hbonds[0][0] == 'O'
        assert hbonds[0][1] == 2.15
        assert hbonds[0][4] == 6
        
    def test3_hbonds_analyzer_compute(self, sample_file3):
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
        
    def test4_hbonds_analyzer_compute(self, sample_file4):
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
        
    def test4_hbonds_analyzer_compute(self, sample_file5):
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
        
    def test5_hbonds_analyzer_compute(self, sample_file6):
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
        
    def test_hbonds_analyzer_render_finalize_output(self, temp_file):
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
        
