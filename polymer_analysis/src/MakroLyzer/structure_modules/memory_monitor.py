"""Opt-in Linux process memory diagnostics, including native-library allocations."""

import csv
import gc
import time
from contextlib import contextmanager


def process_memory_mib():
    """Return current and lifetime peak resident memory from Linux /proc."""
    values = {}
    with open('/proc/self/status') as status:
        for line in status:
            key, _, value = line.partition(':')
            if key in ('VmRSS', 'VmHWM'):
                values[key] = int(value.split()[0]) / 1024
    return values['VmRSS'], values['VmHWM']


class MemoryMonitor:
    def __init__(self, path=None, collect=False):
        self.path = path
        self.collect = collect
        self.file = None

    def __enter__(self):
        if self.path:
            process_memory_mib()  # Fail clearly on unsupported systems before analysis.
            self.file = open(self.path, 'x', newline='')
            self.writer = csv.writer(self.file)
            self.writer.writerow(['frame', 'stage', 'elapsed_seconds', 'rss_mib',
                                  'peak_rss_mib', 'collected_objects'])
            self.started = time.monotonic()
            self.record(-1, 'start')
        return self

    def __exit__(self, *exc):
        if self.file:
            self.file.close()

    def record(self, frame, stage, collected=''):
        if self.file:
            rss, peak = process_memory_mib()
            self.writer.writerow([frame, stage, round(time.monotonic() - self.started, 3),
                                  rss, peak, collected])
            self.file.flush()  # Keep completed samples even if the job is OOM-killed.

    @contextmanager
    def stage(self, frame, name):
        self.record(frame, 'before:' + name)
        try:
            yield
        finally:
            self.record(frame, 'after:' + name)

    def end_frame(self, frame):
        self.record(frame, 'frame_end')
        if self.file and self.collect:
            collected = gc.collect()
            self.record(frame, 'after_gc', collected)
