<!-- Copied/merged from AGENTS.md and README.md. Keep short and actionable. -->
# Guidance for AI coding agents working on lrgsglib

You are an expert Python research engineer. The goal is to produce small, well-tested, readable contributions to lrgsglib (Python + C/C++ extensions).

- Quick build/test workflow (explicit)
  - Create environment: `conda env create -f lrgsgenv.yml && conda activate lrgsgenv`
  - Build native bits and helpers: `make all`
  - Install editable: `pip install -e .`
  - Run tests: `pytest -q` (add/ship tests for any new behavior)

- Project essentials (what to read first)
  - `AGENTS.md` (project-specific agent conventions) — contains the longer persona & dos/don'ts.
  - `README.md` (root) — install and high-level architecture notes.
  - `pyproject.toml` — Python target (>=3.12) and core deps.
  - `Makefile` and `build/` — how C/C++ parts are compiled and configured.

- Architecture & hotspots
  - Languages: Python (primary) with performance-critical C/C++ co-located in `src/lrgsglib/statsys/<Model>/ccore/`.
  - Key Python packages live under `src/lrgsglib/` (utilities, plotting, nx patches, spectral tools).
  - Notebooks live under `ipynb/` and large data under `data/` — avoid committing large binary data.

- Conventions to follow (do not invent these)
  - Prefer small, pure functions, explicit inputs/outputs, and no hidden global state.
  - Determinism: any randomized behavior must accept a `seed` parameter.
  - Use PEP-484 type hints and concise docstrings; include a one-line example when helpful.
  - Always include unit tests (pytest) for new functionality and small doctest examples when useful.
  - Keep diffs minimal and focused; breaking API changes require a migration note.

- C/C++ and performance-critical code
  - Per-model `ccore/` folders contain performance kernels. Do NOT rewrite or remove these files without
    strong tests, benchmarks, and explicit human review.
  - Use the `Makefile` rules and `build/` include files to understand compile flags and paths.

- Integration examples (where to look for patterns)
  - Spectral utilities: `src/lrgsglib/utils/lrg/spectral.py` (examples: computing Laplacian spectra)
  - Plot helpers: `src/lrgsglib/plotlib/` (small utilities for figures used in notebooks)
  - C kernels: `src/lrgsglib/statsys/<Model>/ccore/` (simulators: Ising, Voter, Contact)

- Hard constraints / do nots
  - Do not push/merge to protected branches; always open PRs for human review.
  - Do not add secrets or external network calls.
  - Do not change build, license, or CI defaults without explicit rationale.

- What to include when you change code
  - Short description of intent + files changed.
  - How to reproduce/validate locally (commands used and expected quick checks).
  - If you touch C/C++ code, include a micro-benchmark or unit test showing parity/perf.

If any section is unclear or you want more examples from specific modules, tell me which area (e.g. C build, spectral utils, or plotting conventions) and I will expand this file with targeted examples.
