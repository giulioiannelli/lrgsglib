# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

**lrgsglib** implements the **Laplacian Renormalization Group for Signed Graphs** with Python modules and C/C++ extensions. It provides tools for building signed networks, running renormalization flows, and simulating statistical physics models (Ising, contact process, voter dynamics).

## Build and Installation

**Standalone usage:**
```bash
# Create and activate conda environment
conda env create -f lrgsgenv.yml
conda activate lrgsgenv

# Initialize submodules
git submodule init
git submodule update

# Build C/C++ components and configure environment
make all

# Install in editable mode
pip install -e .
```

**As a submodule** (e.g., in `lrgsglib-ipynb`):
```bash
# From the lrgsglib subdirectory
make all LRGSG_LLIB=$(pwd)/.. CONDA_ENV_NAME=your_env_name
```
- `LRGSG_LLIB` sets the outer project root for data/notebooks/logs
- Without this, paths default to `lrgsglib/` subdirectory
- This configures: data → `outer-project/data`, notebooks → `outer-project/ipynb`, logs → `outer-project/.log`

**Rebuild after changes:**
```bash
make clean
make all
```

## Running Tests

```bash
# From repository root
pytest test/

# Run specific test
pytest test/test_contact_process.py

# Run benchmarks
python test/bench_contact_process_quick.py
python test/bench_mcg_slaplspect_backends.py

# Run diagnostic scripts
python test/diagnose_cp_p0_anomalies.py
```

### Test Folder Organization and Conventions

**CRITICAL: Follow these conventions for all test-related files:**

#### File Naming Conventions
Tests must follow specific naming patterns that clearly indicate their purpose:

- **`test_*.py`** - Unit tests (pytest-compatible)
  - Example: `test_contact_process.py`, `test_lattice2d_spectral.py`
  - For testing specific functionality/modules

- **`bench_*.py`** - Performance benchmarks
  - Example: `bench_contact_process_quick.py`, `bench_mcg_slaplspect_backends.py`
  - Be SPECIFIC: not "benchmark.py" but "bench_<what>_<aspect>.py"
  - Clearly indicate WHAT is being benchmarked in the filename

- **`demo_*.py`** - Demonstration/example scripts
  - Example: `demo_cp2d_relu.py`, `demo_edge_sign_clustering.py`
  - For showcasing features or workflows

- **`diagnose_*.py`** - Diagnostic/debugging scripts
  - Example: `diagnose_cp_p0_anomalies.py`, `diagnose_cp_subcritical_anomalies.py`
  - For investigating specific issues

- **`script_*.py`** - Utility scripts
  - Example: `script_check_border_mode.py`, `script_chl_antiferro_recon.py`
  - For one-off tasks or utilities

#### Directory Structure

```
test/
├── .tmp/                    # ALL temporary/output files go here (gitignored)
│   ├── *.pkl                # Benchmark results
│   ├── *.png                # Generated plots
│   └── *.log                # Log files
├── output/                  # Persistent test outputs (if needed)
├── figures/                 # Test-generated figures for documentation
├── test_*.py                # Unit tests
├── bench_*.py               # Benchmarks
├── demo_*.py                # Demonstrations
└── diagnose_*.py            # Diagnostics
```

#### Output File Management

**IMPORTANT**: Keep test folder clean!

1. **Temporary files** (benchmark results, intermediate data) → `test/.tmp/`
   ```python
   TEST_TMP_DIR = Path(__file__).parent / ".tmp"
   TEST_TMP_DIR.mkdir(exist_ok=True)
   output_file = TEST_TMP_DIR / "my_results.pkl"
   ```

2. **Documentation** (markdown files, reports) → `.agents/devd/<project>/`
   - NEVER put .md files in test folder
   - Agent documentation belongs in `.agents/` tree only

3. **Persistent outputs** (if absolutely necessary) → `test/output/` or `test/figures/`
   - Use sparingly
   - Only for outputs that need to be committed

#### Example: Proper Benchmark Structure

```python
#!/usr/bin/env python
"""
Benchmark MCG_SlaplSpect.py with different backends (scipy vs cupy).

Clear description of what's being tested and why.
Output files stored in test/.tmp/
"""

from pathlib import Path

# Set up output directory
TEST_TMP_DIR = Path(__file__).parent / ".tmp"
TEST_TMP_DIR.mkdir(exist_ok=True)

# ... benchmark code ...

# Save results to .tmp
results_file = TEST_TMP_DIR / "bench_mcg_slaplspect_results.pkl"
```

#### Common Mistakes to Avoid

❌ Vague names: `benchmark.py`, `test.py`, `script.py`
✓ Specific names: `bench_mcg_slaplspect_backends.py`, `test_contact_process.py`

❌ Output files in test root: `test/results.pkl`, `test/output.png`
✓ Output files in .tmp: `test/.tmp/results.pkl`, `test/.tmp/output.png`

❌ Documentation in test folder: `test/BENCHMARK_RESULTS.md`
✓ Documentation in agents folder: `.agents/devd/project/BENCHMARK_RESULTS.md`

❌ Mixed concerns: test + analysis in one file (unless very simple)
✓ Separate files: `bench_X.py` + `bench_X_analyze.py` (if needed)

## Architecture

### Python Library Structure (`src/lrgsglib/`)

- **`config/`** - Configuration, error handling, program argument parsing
  - `progargs/` - Command-line argument definitions for graph types and dynamics (ContactProcess, IsingDynamics, Lattice2D/3D, SignedGraph, etc.)
  - `lrgsg_env.py` - Environment variables from build system

- **`nx_patches/`** - NetworkX extensions for specialized graph types
  - `SignedGraph/`, `Lattice2D/`, `Lattice3D/`, `ErdosRenyi/`, `WeightedGraph/`, etc.
  - `funcs/` - Utility functions (spectral, thresholding, neighbors, lattice operations)

- **`statsys/`** - Python wrappers for statistical physics simulations

- **`plotlib/`** - Plotting utilities specialized for signed graphs and lattices

- **`gt_patches/`** - graph-tool compatibility layer and C++ extensions

- **`utils/`** - Core utilities organized by domain
  - `lrg/` - Laplacian renormalization group tools (spectral analysis, coarse-graining)
  - `basic/` - General utilities
  - `recon/` - Reconstruction algorithms
  - `tools/` - Helper functions

- **`loglib.py`**, **`proglib.py`** - Logging and program execution helpers

### C/C++ Components (`src/lrgsglib/Ccore/`)

Performance-critical code lives here:

- **`statsys/`** - Statistical physics simulators (compiled to `Ccore/bin/`)
  - `RBIsingM/` - Random-bond Ising model (IsingSimulator variants)
  - `voterM/` - Voter model dynamics (VoterSimulator variants)
  - `contactP/` - Contact process (ContactSimulator variants)
  - `signedRw/` - Signed random walks on lattices

- **`SFMT/`** - SIMD-oriented Fast Mersenne Twister RNG

- **`LRGSG_bindynsys.c/h`** - Binary dynamics system utilities

**Build outputs:**
- Binaries: `Ccore/bin/` (e.g., `IsingSimulator0`, `VoterSimulator1`, `ContactSimulator1d`)
- Python extensions: `.so` files via pybind11

### Build System

The Makefile includes modular configuration files from `build/`:
- `lrgsg-paths.mk` - Path definitions (customizable via `LRGSG_LLIB`)
- `cconfig.mk` - C/C++ compiler flags
- `cprogn.mk` - C program build rules
- `conda-config.mk` - Conda environment configuration

Build generates a `.env` file with all paths (`LRGSG_DATA`, `LRGSG_IPYNB`, `LRGSG_LOG`, `LRGSG_CCORE_BIN`, etc.) loaded via `python-dotenv`.

## Common Workflow Patterns

**Running C simulators from Python:**
1. Python scripts (e.g., `src/L2D_IsingDynamics.py`) or kernels (`src/kernels/`) call compiled C binaries
2. Serializers handle data storage (e.g., `L2D_IsingDynamics_Serialiser.py`)
3. Results analyzed in notebooks with plotlib utilities

**Graph creation and analysis:**
1. Create signed graphs using nx_patches (e.g., `Lattice2D.create_graph()`, `ErdosRenyi.signed_erdos_renyi()`)
2. Extract properties (Laplacian spectra, balance ratios) via `utils/lrg/spectral.py`
3. Run dynamics simulations via C programs or Python wrappers
4. Visualize with plotlib

## Code Style Conventions

From `AGENTS.md` and `pyproject.toml`:

- **Python:** PEP 8, PEP 484 type hints, Python ≥3.12
- **Determinism:** Randomized logic must accept a `seed` parameter
- **Modularity:** Prefer pure functions, minimal state, small composable modules
- **Documentation:** Concise docstrings with short runnable examples
- **Testing:** Ship tests with all new code (pytest)
- **Explicit over implicit:** No global state, no wildcard imports

**Formatting tools:**
- `black` (line length 80)
- `isort` (profile: black)
- `mypy` (strict type checking)
- `flake8`

## Key Dependencies

Core: `numpy`, `scipy`, `networkx`, `matplotlib`, `pandas`, `scikit-learn`
Performance: `pybind11`, `cupy` (GPU acceleration)
Graphs: `graph-tool` (from conda)
Specialized: `lmfit`, `powerlaw`, `tqdm`

## Important Notes

- C/C++ code uses SFMT for random number generation
- Build system auto-generates environment configuration (`.env`, `lrgsg_env.py`)
- When working with submodule setup, always rebuild with proper `LRGSG_LLIB` path
- C programs are performance-critical; changes require benchmarks and tests
