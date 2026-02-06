Architecture
============

This document describes the architecture, module organization, and design patterns
used in lrgsglib.

Overview
--------

lrgsglib is organized around four pillars:

1. **Simplicity** - Clear, readable code without unnecessary complexity
2. **Efficiency** - Performance-conscious implementations with appropriate optimizations
3. **Generalizability** - Reusable components that adapt to multiple use cases
4. **Modularity** - Clean separation of concerns with logical submodule organization

Module Organization
-------------------

The library follows a layered architecture where functionality flows from
low-level utilities up through graph operations to high-level programs.

.. code-block:: text

   lrgsglib/
   ├── src/lrgsglib/              ← Core Python library
   │   ├── nx_patches/            ← Graph types (NetworkX extensions)
   │   │   ├── SignedGraph/       ← Base class for all signed graphs
   │   │   ├── Lattice2D/         ← 2D lattice graphs
   │   │   ├── Lattice3D/         ← 3D lattice graphs
   │   │   ├── ErdosRenyi/        ← Random signed graphs
   │   │   ├── MultiplicativeCascade/  ← Hierarchical scale-free graphs
   │   │   ├── DiracLattice/      ← Dirac comb/brush lattices
   │   │   ├── Vicsek/            ← Vicsek fractal graphs
   │   │   └── funcs/             ← Utility functions for graphs
   │   │
   │   ├── statsys/               ← Statistical physics simulations
   │   │   ├── IsingDynamics.py   ← Ising model (Metropolis, SA, PT)
   │   │   ├── ContactProcess.py  ← Contact process (EI, SIR)
   │   │   ├── VoterModel.py      ← Voter model dynamics
   │   │   └── BinDynSys.py       ← Base class for binary dynamics
   │   │
   │   ├── utils/                 ← Utilities organized by domain
   │   │   ├── basic/             ← General utilities (linalg, numeric, I/O)
   │   │   ├── lrg/               ← Domain-specific (spectral, infocomm)
   │   │   └── tools/             ← Data structures (NestedDict, etc.)
   │   │
   │   ├── plotlib/               ← Visualization utilities
   │   ├── config/                ← Configuration and constants
   │   │   └── progargs/          ← Command-line argument definitions
   │   ├── Ccore/                 ← C/C++ performance layer
   │   └── bindings/              ← pybind11 Python-C++ bindings
   │
   ├── src/                       ← Standalone programs
   │   ├── kernels/               ← Reusable computational building blocks
   │   ├── parsers/               ← Argument parsers for programs
   │   ├── L2D_SlaplSpect.py      ← Lattice2D spectral analysis program
   │   ├── MCG_SlaplSpect.py      ← MultiplicativeCascade spectral program
   │   └── ...                    ← Other programs
   │
   └── test/                      ← Tests and benchmarks

Graph Type Hierarchy
--------------------

All graph types inherit from ``SignedGraph``, which extends NetworkX graphs
with signed edge support:

.. code-block:: text

   SignedGraph (base)
   ├── Lattice2D       ← Regular 2D grids (square, triangular, hexagonal)
   ├── Lattice3D       ← Regular 3D grids
   ├── ErdosRenyi      ← Random graphs with signed edges
   ├── MultispectralGraph (abstract base)
   │   ├── MultiplicativeCascadeGraph  ← Scale-free hierarchical
   │   ├── DiracLatticeGraph           ← Dirac comb structures
   │   └── VicsekGraph                 ← Fractal Vicsek structures
   └── WattStrogatz    ← Small-world networks

**Key design principle:** Graph classes store two representations:

- ``gr['G']`` - Graph with integer node labels (for computation)
- ``gr['H']`` - Graph with coordinate node labels (for visualization)

SignedGraph Method Organization
-------------------------------

Large classes are split into multiple files using the **private method pattern**:

.. code-block:: text

   nx_patches/SignedGraph/
   ├── SignedGraph.py     ← Main class definition (~1200 lines)
   ├── _spectral.py       ← Spectral computation methods (~1300 lines)
   ├── _infotheory.py     ← Information theory methods (~200 lines)
   ├── _dynamics.py       ← Dynamics simulation methods (~90 lines)
   ├── _partitioning.py   ← Graph partitioning (~370 lines)
   ├── _topology.py       ← Topological analysis (~200 lines)
   ├── _backend.py        ← Backend management (~620 lines)
   ├── _representations.py ← Matrix representations (~95 lines)
   ├── _loaders.py        ← I/O operations (~170 lines)
   ├── _exports.py        ← Export formats (~220 lines)
   └── _cleaners.py       ← Data cleaning (~50 lines)

Methods in ``_*.py`` files are mixed into the main class at import time.
This pattern keeps individual files manageable while maintaining a unified API.

Program Architecture
--------------------

Programs follow the **kernels + wrappers** pattern:

1. **Kernels** (``src/kernels/``) contain reusable computational building blocks
2. **Parsers** (``src/parsers/``) define command-line arguments
3. **Programs** (``src/*.py``) are thin wrappers combining kernels and parsers

Naming Convention
~~~~~~~~~~~~~~~~~

Programs follow the pattern ``{GRAPH_TYPE}_{COMPUTATION_TYPE}.py``:

- ``L2D_SlaplSpect.py`` - Lattice2D + Signed Laplacian Spectrum
- ``MCG_SlaplSpect.py`` - MultiplicativeCascade + Signed Laplacian Spectrum
- ``L2D_IsingDynamics.py`` - Lattice2D + Ising Model Dynamics

**Serializer programs** add ``_Serialiser.py`` suffix for batch job generation:

- ``L2D_SlaplSpect_Serialiser.py`` - Generate parameter sweep jobs

Minimal Main Program
~~~~~~~~~~~~~~~~~~~~

All main programs follow this template:

.. code-block:: python

   from parsers.L2D_SlaplSpect import parse_arguments, parser
   from kernels.L2D_SlaplSpect import *

   def main():
       args = parse_arguments(parser)
       perform_spectral_calculations(args)  # All logic in kernels
       if args.print_chrono:
           Chronometer.print_all_chronometers()
       if args.verbose:
           print("Done!")

   if __name__ == "__main__":
       main()

This separation enables:

- **Testability** - Import kernels directly in tests
- **Reusability** - Call kernel functions from notebooks
- **Consistency** - All programs have the same structure

Kernel Structure
~~~~~~~~~~~~~~~~

Kernels are organized by:

1. **Graph type** - e.g., ``L2D.py`` (Lattice2D functions)
2. **Computation type** - e.g., ``SlaplSpect.py`` (spectral analysis framework)
3. **Combined** - e.g., ``L2D_SlaplSpect.py`` (combines both)

Example kernel organization:

.. code-block:: python

   # kernels/L2D.py - Low-level lattice operations
   __all__ = [
       "initialize_l2d_dict_args",
       "prepare_lattice",
       "eigV_for_lattice2D",
       "eigv_for_lattice2D",
   ]

   # kernels/SlaplSpect.py - Generic spectral framework
   __all__ = [
       "process_eigen_distribution",
       "save_data",
       "load_existing_data",
       "check_existing_and_needed_averages",
   ]

Backend Strategy
----------------

lrgsglib supports multiple computation backends with automatic selection:

.. list-table:: Backend Options
   :header-rows: 1

   * - Backend
     - Use Case
     - Memory
     - Speed
   * - ``numpy``
     - Small graphs, portability
     - Dense
     - Moderate
   * - ``scipy``
     - Medium graphs, CPU
     - Dense/Sparse
     - Fast
   * - ``cupy``
     - Large graphs, GPU
     - Dense
     - Very Fast

Automatic Backend Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The library implements intelligent fallback:

.. code-block:: python

   def select_optimal_backend(N, requested_backend='cupy'):
       dense_gb = (N ** 2 * 8) / (1024 ** 3)

       if requested_backend == 'cupy' and dense_gb < GPU_VRAM_LIMIT_GB:
           return ('cupy', False, 'denseGPU')
       if dense_gb < CPU_RAM_LIMIT_GB and N < SPARSE_CROSSOVER_N:
           return ('scipy', False, 'denseCPU')
       return ('scipy', True, 'sparseCPU_N-2')

**Fallback chain:** GPU (cupy) → Dense CPU (scipy) → Sparse CPU (scipy)

Python/C Integration
--------------------

Performance-critical code lives in ``Ccore/``:

.. code-block:: text

   Ccore/
   ├── bin/                    ← Compiled binaries
   ├── SFMT/                   ← SIMD Mersenne Twister RNG
   ├── statsys/                ← Simulator source code
   │   ├── RBIsingM/           ← Random-bond Ising model
   │   │   ├── base/           ← Core Ising functions
   │   │   ├── simulatorC/     ← IsingSimulator0, 1a, 1b, 3b, 4b
   │   │   └── storer/         ← Data storage utilities
   │   ├── voterM/             ← Voter model simulators
   │   └── contactP/           ← Contact process simulators
   └── LRGSG_utils.c           ← Common C utilities

**Simulator variants:**

- ``IsingSimulator0`` - Basic Metropolis
- ``IsingSimulator1b`` - Optimized Metropolis with logging
- ``IsingSimulator3b`` - Simulated Annealing
- ``IsingSimulator4b`` - Parallel Tempering

See :doc:`c_extensions` for detailed documentation.

Data Flow
---------

Typical workflow for spectral analysis:

.. code-block:: text

   1. Create Graph
      └─ Lattice2D(side=64, pflip=0.3, geo='sqr')
         └─ Applies sign flips to edges

   2. Compute Spectrum
      └─ sg.compute_laplacian_spectrum_weigV(backend='scipy')
         └─ Builds signed Laplacian: L = D - A_signed
         └─ Computes eigenvalues and eigenvectors

   3. Analyze Results
      └─ utils/lrg/infocomm.py functions
         └─ compute_entropy_observables_from_eigenvalues()

   4. Visualize
      └─ plotlib utilities
         └─ plot_eigenvalue_distribution()

Typical workflow for dynamics simulation:

.. code-block:: text

   1. Create Graph
      └─ Lattice2D(side=64, pflip=0.3)

   2. Initialize Dynamics
      └─ IsingDynamics(sg=lattice, T=2.0, runlang='C1b')
         └─ init_ising_dynamics()

   3. Run Simulation (C backend)
      └─ Export graph data to binary files
      └─ Call IsingSimulator1b binary
      └─ Read results from stdout

   4. Analyze
      └─ Access ising.magn, ising.ene, ising.s

Configuration System
--------------------

Configuration is centralized in ``config/``:

.. code-block:: text

   config/
   ├── __init__.py       ← Imports from const.py, errwar.py, funcs.py
   ├── const.py          ← Constants (paths, file extensions)
   ├── errwar.py         ← Custom exceptions and warnings
   ├── funcs.py          ← Configuration functions
   ├── lrgsg_env.py      ← Auto-generated environment variables
   └── progargs/         ← Program argument definitions
       ├── ContactProcess.py
       ├── IsingDynamics.py
       ├── Lattice2D.py
       └── ...

Environment variables (from ``.env`` file):

- ``LRGSG_ROOT`` - Library root directory
- ``LRGSG_DATA`` - Data storage directory
- ``LRGSG_IPYNB`` - Notebooks directory
- ``LRGSG_LOG`` - Log files directory
- ``LRGSG_CCORE_BIN`` - C binary directory

Design Patterns Summary
-----------------------

.. list-table::
   :header-rows: 1

   * - Pattern
     - Where Used
     - Purpose
   * - Private Method Files
     - ``nx_patches/*/``
     - Split large classes into manageable files
   * - Kernels + Wrappers
     - ``src/``, ``src/kernels/``
     - Separate computation from CLI
   * - Callback-based Framework
     - ``SlaplSpect.py``
     - Generic framework with graph-specific callbacks
   * - Backend Fallback
     - ``_spectral.py``, ``_backend.py``
     - Graceful degradation across backends
   * - Checkpoint/Resume
     - ``SlaplSpect.py``
     - Resume long computations

Function Placement Guidelines
-----------------------------

**When to place in ``kernels/``:**

- Reusable computational building blocks
- Graph-specific operations (e.g., lattice initialization)
- Generic algorithms that work across graph types

**When to place in ``utils/``:**

- General-purpose utilities (math, I/O)
- Information-theoretic calculations
- Data structure helpers

**When to place in ``nx_patches/``:**

- Graph class definitions
- Methods that belong to graph objects
- Internal implementation files (``_*.py``)

See Also
--------

- :doc:`build_system` - Build process details
- :doc:`c_extensions` - C/C++ integration
- :doc:`style_guide` - Code conventions
