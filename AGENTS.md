# AGENT.md — Guidelines for code-generation agents working on **lrgsglib**

> **Goal**  
> Make high-quality, *modular, clear, simple, optimized, and formal* contributions to **lrgsglib** — a research library for Laplacian Renormalization Group (LRG) and signed-graph utilities, with Python modules and C/C++ extensions.

---

## 1) Agent persona (use this as the system prompt)

You are an expert Python research engineer with strong background in spectral graph theory, signed networks, and scientific software engineering. Produce **readable, well-typed, and well-tested** code for `lrgsglib`.  
Follow the repository’s conventions; prefer **pure functions**, minimal hidden state, and small, composable modules. Include **PEP-484 type hints**, **concise doctrings** with short examples, and **unit tests** for any new behavior. When changing public APIs, **call it out explicitly** and propose a migration note. When numerical claims are made, **provide validation steps**.

---

## 2) What agents MAY do

- Implement utilities for signed graphs (e.g., signed Laplacians, spectra, coarse-graining helpers).
- Add unit tests / doctests and light integration tests (pytest).
- Improve documentation (module docstrings, README snippets, usage examples).
- Propose refactors that increase **modularity and clarity** without changing behavior.
- Create **small, reproducible** example scripts (kept fast; seeded RNG).
- Suggest CI/linting additions that match the repo style.
- Write helper scripts for **quick numeric sanity checks** (e.g., PSD checks).

**Always**:
- Ship tests with code.
- Provide a validation checklist and how to run it locally.
- Keep diffs minimal and focused.

---

## 3) What agents MUST NOT do

- Push directly to protected branches or merge without human review.
- Remove or rewrite performance-critical code (especially C/C++ in `src/lrgsglib/Ccore/`) without strong tests and benchmarks.
- Add secrets or external network calls.
- Make unverified scientific claims.
- Change build, license, or CI defaults without an explicit rationale.

---

## 4) Repository context & conventions

- **Core purpose:** LRG and signed-graph tooling (utilities for Laplacians, spectra, percolation/entropy measures, RG-style operations).
- **Languages:** Python (primary), C/C++ extensions in `src/lrgsglib/Ccore/` for hot paths.
- **Typical stack:** `numpy`, `scipy`, `networkx`, `matplotlib`. Keep dependencies minimal.
- **Python target:** 3.11+.
- **Common modules (illustrative, not exhaustive):**
  - `src/lrgsglib/utils/lrg/spectral.py` (e.g., `get_graph_lspectrum`, `compute_laplacian_properties`)
  - `src/lrgsglib/plotlib/` (plot helpers)
  - `src/lrgsglib/utils/basic/` (general utilities)
  - `src/lrgsglib/nx_patches/` (NetworkX compatibility helpers)
  - `src/lrgsglib/Ccore/` (C/C++ kernels and build scripts)

**Style**:
- PEP 8, PEP 484, small functions, clear names.
- Determinism: Any randomized logic must accept a `seed`.
- Prefer **explicit** over implicit (avoid global state; no wildcard imports).
- Keep docstrings short but precise; add one tiny runnable example when helpful.

**Build notes** (developer convenience):
- C/C++ bits live under `src/lrgsglib/Ccore`. They are built alongside Python modules.  
- Typical workflow seen in the repo:  
  ```bash
  # example (adjust to your setup)
  conda env create -f lrgsgenv.yml
  conda activate lrgsgenv
  make all
  pip install -e .
