# Investigating memory growth

For a structure-analysis run on Linux, append `--memory-log memory.csv` to
the usual MakroLyzer command. Use a new filename for each run; existing logs
are never overwritten. No extra dependency is required.

The CSV records each analyzed frame after loading, after graph preparation,
before and after each analyzer, and at frame end. Rows are flushed immediately
so completed samples remain available if the scheduler kills the process.
Frames are numbered from zero, like analysis outputs.

- `rss_mib`: current resident process memory in MiB, including native libraries.
- `peak_rss_mib`: highest resident memory reached since process startup. This
  value never decreases, even when temporary memory has been released.
- `collected_objects`: objects reclaimed during optional forced collection.

Compare `frame_end` RSS across frames after initial caches have warmed up.
Persistent growth warrants investigation, but does not by itself prove a leak:
allocators may retain freed memory for reuse. Compare before/after rows for
each analyzer to locate growth. A peak increase with little change in current
RSS suggests a temporary allocation. Measurements at stage boundaries cannot
describe every intermediate allocation; the lifetime peak helps expose spikes.
The last `before:` row without a corresponding `after:` row can identify the
analyzer active when the job was killed.

Repeat with `--memory-log memory_gc.csv --memory-gc` to force Python garbage
collection after each analyzed frame. If this stops the growth, delayed cleanup
of unreachable Python objects is a candidate. Collection may not lower RSS,
and cannot free live objects or native allocations still held by a library.
This option changes cleanup timing and may slow analysis; use it for comparison.

Also compare runs with surface analysis (`-sA`) and convex-hull solvent analysis
(`-ConvHullSol ...`) separately. Keep the trajectory, other options, and scheduler
memory allowance consistent. Sampling fewer frames with `-nth` reduces runtime,
but can miss a geometry-dependent spike and changes the number of allocations.

The monitor measures this process, not total memory used by all tasks in a
scheduler job. Synthetic tests verify logging and cycle reclamation; confirming
the real cause requires running the diagnostic on the affected trajectory.
