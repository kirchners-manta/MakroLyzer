import csv
import gc
import weakref

import pytest

from MakroLyzer.structure_modules.memory_monitor import MemoryMonitor, process_memory_mib


def test_current_and_peak_memory():
    rss, peak = process_memory_mib()
    assert 0 < rss <= peak


def test_log_survives_analyzer_exception(tmp_path):
    path = tmp_path / 'memory.csv'
    with pytest.raises(RuntimeError):
        with MemoryMonitor(path) as monitor:
            with monitor.stage(5, 'FailingAnalyzer'):
                raise RuntimeError('analysis failed')
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    assert [row['stage'] for row in rows] == [
        'start', 'before:FailingAnalyzer', 'after:FailingAnalyzer']
    assert rows[-1]['frame'] == '5'
    # A diagnostic rerun must not silently overwrite previous evidence.
    with pytest.raises(FileExistsError):
        with MemoryMonitor(path):
            pass


def test_forced_collection_reclaims_unreachable_cycles(tmp_path):
    class Cycle:
        pass

    was_enabled = gc.isenabled()
    gc.disable()
    try:
        with MemoryMonitor(tmp_path / 'gc.csv', collect=True) as monitor:
            obj = Cycle()
            obj.self = obj
            reference = weakref.ref(obj)
            del obj
            assert reference() is not None
            monitor.end_frame(0)
            assert reference() is None
    finally:
        if was_enabled:
            gc.enable()


def test_disabled_monitor_does_not_force_collection(monkeypatch):
    def unexpected():
        raise AssertionError('normal analysis must not force collection')

    monkeypatch.setattr(gc, 'collect', unexpected)
    with MemoryMonitor(collect=True) as monitor:
        with monitor.stage(0, 'Analyzer'):
            pass
        monitor.end_frame(0)
