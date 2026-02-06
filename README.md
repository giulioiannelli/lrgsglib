# LRGSGlib

LRGSGlib is a set of Python modules and C/C++ extensions implementing the theoretical tools of the **Laplacian Renormalisation Group for Signed Graphs**. It provides utilities for building signed networks, running renormalisation flows and simulating statistical physics models such as the Ising, contact process and voter dynamics. Additional helpers for plotting, logging and networking patches are also included.

The C sources for the performance critical parts live in `src/lrgsglib/Ccore` and can be built along with the Python modules.

## Installation

1. **Clone the repository and initialise submodules**
   ```bash
   git clone https://github.com/giulioiannelli/lrgsglib
   cd lrgsglib
   git submodule init
   git submodule update
   ```

2. **Create the conda environment** (default name: `lrgsgenv`)
   ```bash
   conda env create -f lrgsgenv.yml
   conda activate lrgsgenv
   ```

3. **Build the project**

   **Standalone usage:**
   ```bash
   make all
   ```

   **Submodule usage** (e.g., in `lrgsglib-ipynb`):
   Configure paths relative to the outer project instead of the `lrgsglib` subdirectory:
   ```bash
   # From the lrgsglib subdirectory
   make all LRGSG_LLIB=$(pwd)/.. CONDA_ENV_NAME=your_env_name
   ```

   **What this does:**
   - Sets `LRGSG_LLIB` to the outer project root
   - Configures data folder as `outer-project/data` (instead of `lrgsglib/data`)
   - Configures notebooks as `outer-project/ipynb` (instead of `lrgsglib/ipynb`)
   - Configures logs as `outer-project/.log` (instead of `lrgsglib/.log`)
   - Keeps all library source code paths relative to `lrgsglib/`
   - Uses your custom conda environment name instead of the default `lrgsgenv`

   **Example for lrgsglib-ipynb:**
   ```bash
   cd lrgsglib
   make all LRGSG_LLIB=$(pwd)/.. CONDA_ENV_NAME=lrgsgnb
   ```

4. **Install in editable mode**
   ```bash
   pip install -e .
   ```

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

## Running Tests

```bash
cd lrgsglib
pytest test/
```

## More Information

- **C programs and models**: See `src/lrgsglib/Ccore/README.md`
- **Agent workspace**: See `.agents/00-START-HERE.md` for development context
- **Examples**: See `docs/source/user_guide/examples.rst` for complete workflows
