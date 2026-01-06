import pytest
import tempfile
from pathlib import Path
import numpy as np

from unittest.mock import Mock, patch

from MakroLyzer.input_handling import readXYZ
from MakroLyzer import graphs

from MakroLyzer.modify_modules.structureModifierBase import ModifyOutputHandler, StructureModifier
from MakroLyzer.modify_modules.Patterns import PatternModifier

# OutputHandler Tests # ------------------------------------------------------------------
class TestModifyOutputHandler:
    """Test the ModifyOutputHandler class."""
    
    # SetUp fixtures -----------------------------------------
    
    @pytest.fixture
    def temp_file(self):
        """Fixture that creates a temporary file path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir) / "output.csv"
    
    # ---------------------------------------------------------
    
    def test_TestModifyOutputHandler_init(self, temp_file):
        """Test initialization of ModifyOutputHandler."""
        handler = ModifyOutputHandler(temp_file)
        assert handler.file_path == temp_file
        
    # ---------------------------------------------------------
    
    def test_TestModifyOutputHandler_write_output(self, temp_file):
        """Test writing a file."""
        handler = ModifyOutputHandler(temp_file)
        header = "21\nx,y,z,ananas"
        
        class _Nodes:
            def __init__(self, mapping, order):
                self._mapping = mapping
                self._order = order

            def __call__(self):
                return list(self._order)

            def __getitem__(self, key):
                return self._mapping[key]

        mock_node1 = object()
        mock_node2 = object()
        mapping = {
            mock_node1: {'index': 0, 'element': 'C', 'x': 1.0, 'y': 2.0, 'z': 3.0, 'ananas': 'yes'},
            mock_node2: {'index': 1, 'element': 'O', 'x': 4.0, 'y': 5.0, 'z': 6.0, 'ananas': 'no'}
        }
        mock_graph = type('G', (), {})()
        mock_graph.nodes = _Nodes(mapping, [mock_node1, mock_node2])
        attributes = ['ananas']
        handler.write_output(header, mock_graph, attributes)
        
        assert temp_file.exists()
        content = temp_file.read_text()
        assert content == "21\nx,y,z,ananas\nC         1.0000    2.0000    3.0000    yes\nO         4.0000    5.0000    6.0000    no\n"
        
# StructureModifier Tests # ------------------------------------------------------------------

class DummyModifier(StructureModifier):
    """A dummy subclass of StructureModifier for testing purposes."""
    
    def compute(self, graph):
        mock_graph = Mock()
        return mock_graph
    
    def render_output(self, graph, frame_idx):
        pass
    
class TestStructureModifier:
    """Test the StructureModifier base class and its subclass."""
    
    def test_StructureModifier_init(self):
        """Test initialization of StructureModifier."""
        handler = Mock()
        modifier = DummyModifier(output_handler=handler)
        assert modifier.output_handler == handler
        assert modifier.frame_number == 0
        
    def test_StructureModifier_compute_and_render(self):
        """Test compute and render_output methods of DummyModifier."""
        handler = Mock()
        modifier = DummyModifier(output_handler=handler)
        
        mock_graph = Mock()
        result = modifier.compute(mock_graph)
        assert isinstance(result, Mock)
        
        modifier.render_output(result, frame_idx=0)
        
# PatternModifier Tests # ------------------------------------------------------------------

class TestPatternModifier:
    """Test the PatternModifier class."""
    
    # SetUp fixtures -----------------------------------------
    
    @pytest.fixture
    def temp_file(self):
        """Fixture that creates a temporary file path."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir) / "output.csv"
            
    @pytest.fixture
    def sample_data1(self):
        return 'test_structures/Patterns/01.xyz'

    @pytest.fixture
    def pattern_data1(self):
        return 'test_structures/Patterns/01pattern.txt'
    
    # ---------------------------------------------------------
    
    def test_PatternModifier_init(self, pattern_data1):
        """Test initialization of PatternModifier."""
        handler = Mock()
        modifier = PatternModifier(pattern_data1, output_handler=handler)
        assert modifier.output_handler == handler

    def test_PatternModifier_read_patterns(self, temp_file):
        """_read_patterns should parse a valid pattern file into a DataFrame."""
        pattern_file = temp_file.parent / "pat.txt"
        # first line: list of lists, second line: element
        pattern_file.write_text("[[1,2],[3,4]]\nC\n")

        modifier = PatternModifier(str(pattern_file))
        df = modifier.patterns
        assert 'pattern' in df.columns
        assert 'element' in df.columns
        assert df['pattern'].values[0] == [[1, 2], [3, 4]]
        assert df['element'].values[0] == 'C'

    def test_PatternModifier_compute(self):
        """compute should call graph.find_and_tag_patterns with patterns and element."""
        # prepare a small pattern file
        import tempfile as _tf
        p = _tf.NamedTemporaryFile(mode='w', delete=False)
        p.write('[[10,20]]\nO\n')
        p.close()

        modifier = PatternModifier(p.name)

        class FakeGraph:
            def __init__(self):
                self.called = False
                self.args = None

            def find_and_tag_patterns(self, pattern, element):
                self.called = True
                self.args = (pattern, element)

        g = FakeGraph()
        res = modifier.compute(g)
        assert res is g
        assert g.called is True
        assert g.args[0] == [[10, 20]]
        assert g.args[1] == 'O'

    def test_PatternModifier_render_output(self, temp_file):
        """render_output should call output_handler.write_output with header and fragmentID attribute."""
        # pattern file
        pattern_file = temp_file.parent / "pat2.txt"
        pattern_file.write_text("[[1]]\n\n")

        mock_handler = Mock()
        modifier = PatternModifier(str(pattern_file), output_handler=mock_handler)

        class G:
            def number_of_nodes(self):
                return 2

        g = G()
        modifier.render_output(g, frame_idx=5)

        expected_header = f"2\nelement x          y         z          FragmentID"
        mock_handler.write_output.assert_called_once_with(expected_header, g, attributes=['fragmentID'])

# AminoSaturationModifier Tests # ------------------------------------------------------------------

class TestAminoSaturationModifier:
    """Tests for the AminoSaturationModifier."""
    
    # SetUp fixtures -----------------------------------------

    @pytest.fixture
    def temp_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir) / "output.csv"
            
    # ---------------------------------------------------------

    def test_AminoSaturationModifier_init(self):
        handler = Mock()
        from MakroLyzer.modify_modules.AminoSaturation import AminoSaturationModifier
        mod = AminoSaturationModifier(output_handler=handler)
        assert mod.output_handler == handler

    def test_AminoSaturationModifier_compute(self):
        from MakroLyzer.modify_modules.AminoSaturation import AminoSaturationModifier

        class G:
            def __init__(self):
                self.saturated = False

            def saturate(self):
                self.saturated = True

        g = G()
        mod = AminoSaturationModifier()
        res = mod.compute(g)
        assert res is g
        assert g.saturated is True

    def test_AminoSaturationModifier_render_output(self, temp_file):
        from MakroLyzer.modify_modules.AminoSaturation import AminoSaturationModifier
        mock_handler = Mock()
        mod = AminoSaturationModifier(output_handler=mock_handler)

        class G:
            def number_of_nodes(self):
                return 4

        g = G()
        mod.render_output(g, frame_idx=0)
        expected = f"4\nelement  x         y        z"
        mock_handler.write_output.assert_called_once_with(expected, g)
        
    
# SubgraphPrintModifier Tests # ------------------------------------------------------------------

class TestSubgraphPrintModifier:
    """Tests for the SubgraphPrintModifier."""

    def test_SubgraphPrintModifier_init(self):
        from MakroLyzer.modify_modules.SubgraphPrint import SubgraphPrintModifier
        handler = Mock()
        mod = SubgraphPrintModifier(output_handler=handler)
        assert mod.output_handler == handler
        assert mod.frame_number == 0

    def test_SubgraphPrintModifier_compute(self):
        from MakroLyzer.modify_modules.SubgraphPrint import SubgraphPrintModifier

        class G:
            def get_subgraphs(self):
                g1 = Mock()
                g2 = Mock()
                return [g1, g2]

        mod = SubgraphPrintModifier()
        res = mod.compute(G())
        assert isinstance(res, list)
        assert len(res) == 2

    def test_SubgraphPrintModifier_render_output(self):
        from MakroLyzer.modify_modules.SubgraphPrint import SubgraphPrintModifier
        mock_handler = Mock()
        mod = SubgraphPrintModifier(output_handler=mock_handler)

        s1 = Mock()
        s1.number_of_nodes.return_value = 2
        s2 = Mock()
        s2.number_of_nodes.return_value = 3
        subgraphs = [s1, s2]

        mod.render_output(subgraphs, frame_idx=0)

        expected1 = f"2\nelement  x         y        z"
        expected2 = f"3\nelement  x         y        z"
        mock_handler.write_output.assert_any_call(expected1, s1, ID=1)
        mock_handler.write_output.assert_any_call(expected2, s2, ID=2)
        assert mock_handler.write_output.call_count == 2

# FunctionalizePEModifier Tests # ------------------------------------------------------------------
class TestFunctionalizePEModifier:
    """Tests for the FunctionalizePEModifier."""
    
    def test_FunctionalizePEModifier_init(self):
        from MakroLyzer.modify_modules.FunctionalizePE import FunctionalizePEModifier
        handler = Mock()
        mod = FunctionalizePEModifier(15, 'CO', output_handler=handler)
        assert mod.output_handler == handler
        assert mod.frame_number == 0
        assert mod.percentage == 15
        assert mod.func_type == 'CO'
        
    def test_FunctionalizePEModifier_compute(self):
        from MakroLyzer.modify_modules.FunctionalizePE import FunctionalizePEModifier
        mock_graph = Mock()
        mod = FunctionalizePEModifier(25, 'CO')
        res = mod.compute(mock_graph)
        assert res is mock_graph
        mock_graph.functionalizePE.assert_called_once_with(25, 'CO')

    def test_FunctionalizePEModifier_render_output(self):
        from MakroLyzer.modify_modules.FunctionalizePE import FunctionalizePEModifier
        mock_handler = Mock()
        mod = FunctionalizePEModifier(10, 'CO', output_handler=mock_handler)

        class G:
            def number_of_nodes(self):
                return 5

        g = G()
        mod.render_output(g, frame_idx=0)
        expected = f"5\nelement  x         y        z"
        mock_handler.write_output.assert_called_once_with(expected, g)
