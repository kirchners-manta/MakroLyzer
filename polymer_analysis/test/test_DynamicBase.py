"""Exercise per-frame rendering and separate final computation/output."""
from unittest.mock import Mock

from MakroLyzer.dynamic_modules.dynamicBase import DynamicAnalyzer


class PerFrameAnalyzer(DynamicAnalyzer):
    def compute(self, graph):
        return graph

    def render_output(self, data, frame_idx):
        self.output_handler.append_row(f'{frame_idx},{data}')


class AccumulatingAnalyzer(PerFrameAnalyzer):
    def compute(self, graph):
        return None

    def finalize(self):
        self.results = {'samples': self.frame_number}
        return self.results

    def finalize_output(self, header=None):
        if self.output_handler is not None:
            self.output_handler.write_csv('samples', [str(self.results['samples'])])


def test_frame_index_reaches_output_and_zero_is_rendered():
    handler = Mock(mode='collect')
    analyzer = PerFrameAnalyzer(handler)
    assert analyzer.run(0., frame_idx=12) == 0.
    handler.append_row.assert_called_once_with('12,0.0')
    assert analyzer.frame_number == 1
    assert analyzer.finalize() is None
    handler.finalize.assert_not_called()
    analyzer.finalize_output('frame,value')
    handler.finalize.assert_called_once_with('frame,value')


def test_accumulated_calculation_and_output_are_separate():
    handler = Mock(mode='collect')
    analyzer = AccumulatingAnalyzer(handler)
    assert analyzer.run(None, frame_idx=0) is None
    assert analyzer.run(None, frame_idx=5) is None
    assert handler.method_calls == []
    assert analyzer.finalize() == {'samples': 2}
    assert handler.method_calls == []
    analyzer.finalize_output()
    handler.write_csv.assert_called_once_with('samples', ['2'])
    handler.append_row.assert_not_called()


def test_results_available_without_output_handler():
    analyzer = PerFrameAnalyzer()
    assert analyzer.run(2., 3) == 2.
    assert analyzer.finalize() is None
    analyzer.finalize_output()
    accumulator = AccumulatingAnalyzer()
    accumulator.run(None, 0)
    assert accumulator.finalize() == {'samples': 1}
    accumulator.finalize_output()


def test_streaming_output_is_not_flushed_as_collected_rows():
    handler = Mock(mode='streaming')
    analyzer = PerFrameAnalyzer(handler)
    analyzer.run(1., 5)
    analyzer.finalize()
    analyzer.finalize_output()
    handler.append_row.assert_called_once_with('5,1.0')
    handler.finalize.assert_not_called()
