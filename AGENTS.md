# MakroLyzer

## Session summaries

- At the beginning of each new session, read `SESSION_SUMMARY.md` before
  starting project work.
- Before the final response after substantive work or planning, refresh
  `SESSION_SUMMARY.md` so a handoff exists if the session closes. Keep one
  concise current summary, normally within 500 words, with the date,
  agreed goals, completed changes and verification, unresolved decisions,
  and the next task. Distinguish proposals from implemented behavior and
  replace superseded conclusions with their corrected versions.
- Include linked references for literature mentioned by the user or the
  assistant, with a short indication of relevance. Carry forward references
  and unresolved plans needed for the next session; link longer notes for
  details.

## Development

- Install from `polymer_analysis/`, following the README:
  `python -m pip install -e ".[test]"`.
  Include the `builder` extra when working on builders requiring RDKit.
- Always run tests with `polymer_analysis/` as the working directory:
  `python -m pytest test`.
  Test fixtures use relative `test_structures/` paths that require this directory.
  Select the relevant test files for focused changes.
- Test realistic pipeline behavior using inputs and intermediate states
  reachable through MakroLyzer's actual entry points and processing steps.
  Cover normal workflows and reachable failures, including invalid user
  inputs at the validation boundary. Avoid contrived internal states or
  cases that earlier pipeline stages already prevent; mocks should preserve
  the real pipeline's contracts.
- Prefer the smallest complete change that solves the requested problem.
  Reuse existing project functionality and dependencies. Add abstractions,
  options, or dependencies only when the task needs them. Keep scientific
  calculations readable and preserve numerical correctness and necessary
  validation.
- Always accompany new code with clear, plain-language comments that help
  the maintainer understand its purpose and how it works. Explain the main
  processing steps, non-obvious choices, and scientific assumptions; include
  units where relevant. Use docstrings to explain new functions and classes,
  including their inputs and outputs. Focus on useful explanations rather
  than restating each line of Python.
- Match the surrounding Python style. The existing pre-commit
  configuration is `polymer_analysis/pre-commit-config.yaml`.

## Architecture

- `polymer_analysis/src/MakroLyzer/MakroLyzer.py` dispatches CLI work
  to structure analysis, dynamic analysis, modification, and building.
- `input_handling/` handles CLI arguments and structure/trajectory input.
- `graphs.py` provides `GraphManager`, a NetworkX graph subclass.
- Extend analyzers, modifiers, and builders through their existing
  base classes and registries in the corresponding module directories.
- User documentation is Sphinx/reStructuredText under `docs/source/`.

## Agent skills

### Issue tracker

Use local Markdown under `.scratch/`. Before creating specs or tickets,
read `docs/agents/issue-tracker.md`.

### Triage labels

Use the default status vocabulary. Before triaging tickets,
read `docs/agents/triage-labels.md`.

### Domain docs

Use the single-context layout. Before exploring domain concepts or
proposing design changes, read `docs/agents/domain.md`.
