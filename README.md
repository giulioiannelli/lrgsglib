# LRGSGlib

LRGSGlib is a set of Python modules and C/C++ extensions implementing the theoretical tools of the **Laplacian Renormalisation Group for Signed Graphs**. It provides utilities for building signed networks, running renormalisation flows and simulating statistical physics models such as the Ising, contact process and voter dynamics. Additional helpers for plotting, logging and networking patches are also included.

The C sources for the performance critical parts live in `src/lrgsglib/Ccore` and can be built along with the Python modules.

## Installation

Choose an installation path depending on your needs:

| Path | C Backends | graph-tool | Compiler needed | Setup |
|------|-----------|------------|-----------------|-------|
| **pip install** | No (Python/numba fallback) | No | No | `pip install .` |
| **pip wheel** (pre-compiled) | Yes (bundled) | No | No | `pip install lrgsglib` |
| **pixi** (recommended) | Yes (from source) | Optional | Provided by pixi | `pixi install` |
| **conda** | Yes (from source) | Optional | Provided by conda | `conda env create` |

### Option 1: pip install (lightweight, no compiler)

```bash
git clone --recursive https://github.com/giulioiannelli/lrgsglib
cd lrgsglib
pip install .
```

All spectral analysis, entropy, graph construction and plotting work. Dynamics
simulations run via Python/numba backends (`runlang="Python"`). If C binaries
are not found, the library warns and falls back to Python automatically.

### Option 2: pixi (recommended for development)

```bash
git clone --recursive https://github.com/giulioiannelli/lrgsglib
cd lrgsglib
pixi install              # Creates environment with GCC, pybind11, etc.
pixi run build            # Compiles C code + installs editable
pixi run test-quick       # Verify everything works
```

For graph-tool support: `pixi install -e full`

### Option 3: conda (still supported)

```bash
git clone --recursive https://github.com/giulioiannelli/lrgsglib
cd lrgsglib
conda env create -f lrgsgenv.yml
conda activate lrgsgnb
pip install -e . --no-build-isolation
```

### Option 4: Makefile (manual build)

For full control over the build process:

```bash
conda activate lrgsgnb
make all                  # Standalone
# or
make all LRGSG_LLIB=$(pwd)/.. CONDA_ENV_NAME=lrgsgnb  # As submodule
pip install -e .
```

### Docker (for reproducibility)

A Dockerfile is provided for reproducible environments, CI, or sharing with
collaborators who don't want to install anything locally:

```bash
docker build -t lrgsglib .
docker run -it lrgsglib python -c "import lrgsglib; print('OK')"
docker run -v ./scripts:/scripts lrgsglib python /scripts/my_analysis.py
```

This is not a development setup — use pixi or conda for active work.

## Documentation

Full documentation is available in `docs/` and can be built with Sphinx:

```bash
cd docs
make html
# Open docs/build/html/index.html in your browser
```

The documentation includes:

- **User Guide** - Getting started, spectral analysis, dynamics simulations, and complete examples
- **API Reference** - Full reference for all modules including:
  - `lrgsglib.graphs` - Unified graph interface with multi-engine support (NetworkX, graph-tool)
  - `lrgsglib.statsys` - Statistical physics models (Ising, Contact Process, Voter)
  - `lrgsglib.utils.lrg` - Spectral analysis and entropy computation
  - `lrgsglib.plotlib` - Visualization utilities
- **Theory** - Mathematical background on signed graphs and the Laplacian RG
- **Developer Guide** - Contributing guidelines and architecture overview

### Key Features

- **Multi-engine graph support**: Create graphs with NetworkX or graph-tool backends
  ```python
  from lrgsglib.graphs import Lattice2D, ErdosRenyi, GraphOfGraphs

  # Explicit engine selection
  lat_nx = Lattice2D(side1=100, geo='sqr', engine='nx')
  lat_gt = Lattice2D(side1=100, geo='sqr', engine='gt')

  # Hierarchical graph structures
  gog = GraphOfGraphs(
      base_graph_type='Lattice2D', base_params={'side1': 10},
      fiber_graph_type='ErdosRenyi', fiber_params={'n': 20, 'p': 0.2}
  )
  ```

- **Spectral entropy analysis**: Compute von Neumann entropy and spectral dimension
  ```python
  from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
  from lrgsglib.utils.lrg.infocomm import compute_entropy_observables_from_eigenvalues

  L, eigenvalues = get_graph_lspectrum(graph.gr['G'])
  entropy, heat, var, tau = compute_entropy_observables_from_eigenvalues(eigenvalues, N)
  ```

- **High-performance C backends**: ~100x speedup for dynamics simulations
  ```python
  from lrgsglib.statsys.IsingDynamics import IsingDynamics

  ising = IsingDynamics(sg=lattice, T=2.0, steps=10000, runlang='C1b')
  ```

- **Graceful degradation**: If C binaries are missing, dynamics automatically fall back to Python with a warning. No manual configuration needed.

## Running Tests

```bash
cd lrgsglib
pytest test/
pytest test/ --quick   # Fast subset (~5s)
```

## More Information

- **C programs and models**: See `src/lrgsglib/Ccore/README.md`
- **Agent workspace**: See `.agents/00-START-HERE.md` for development context
- **Examples**: See `docs/source/user_guide/examples.rst` for complete workflows
