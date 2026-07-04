Installation
============

lrgsglib can be installed either as a standalone library or as a git submodule within a larger project. The library requires Python 3.12+ and includes C/C++ extensions that need to be compiled.

Requirements
------------

**System Requirements:**

- Python 3.12 or later
- C/C++ compiler (GCC or Clang)
- conda (recommended) or Python virtual environment
- Git (for cloning and submodule management)

**Python Dependencies:**

Core dependencies are automatically installed via pip:

- numpy >= 1.19.0
- scipy >= 1.7.0
- networkx >= 2.6.0
- matplotlib >= 3.5.0
- pandas >= 1.3.0
- scikit-learn >= 1.0.0
- tqdm >= 4.60.0
- lmfit >= 1.0.0
- powerlaw >= 1.4.0
- joblib >= 1.0.0
- numba >= 0.58.0
- python-dotenv >= 1.0.0
- plotly >= 5.0.0
- requests >= 2.25.0

Optional dependencies:

- cupy >= 9.0.0 (for GPU acceleration)
- graph-tool (from conda, for high-performance graph operations)

Build-time requirements (installed automatically by ``pip`` in an isolated
build environment; not runtime dependencies):

- scikit-build-core >= 0.11
- pybind11 >= 2.12
- setuptools_scm >= 8

Standalone Installation
------------------------

For using lrgsglib as an independent library:

1. **Clone the repository and initialize submodules:**

   .. code-block:: bash

      git clone https://github.com/giulioiannelli/lrgsglib
      cd lrgsglib
      git submodule init
      git submodule update

2. **Create and activate the conda environment:**

   The default environment name is ``lrgsgenv``:

   .. code-block:: bash

      conda env create -f lrgsgenv.yml
      conda activate lrgsgenv

3. **Install in editable mode:**

   ``pip`` uses the ``scikit-build-core`` backend, so this single step
   compiles the C/C++ (CMake + pybind11) extensions automatically — no
   separate build step is required:

   .. code-block:: bash

      pip install -e .

   For development with additional tools:

   .. code-block:: bash

      pip install -e ".[dev,jupyter,docs]"

.. note::

   A legacy ``make`` route exists for maintainers who want to compile the
   C/C++ simulators outside of ``pip`` (``make all`` in the repository root).
   It is **not** required to install or use the library. In particular, do
   **not** run ``make env-config``: it overwrites
   ``src/lrgsglib/config/lrgsg_env.py`` with machine-pinned absolute paths.
   To relocate output directories, set the ``LRGSG_DATA`` (and optionally
   ``LRGSG_LLIB``/``LRGSG_LOG``/``LRGSG_IPYNB``) environment variables — or a
   ``.env`` file — instead. Data and log directories are created lazily on
   first write.

Submodule Installation
-----------------------

For using lrgsglib as a git submodule within a larger project (e.g., ``lrgsglib-ipynb``):

1. **Add as a submodule** (if not already present):

   .. code-block:: bash

      git submodule add https://github.com/giulioiannelli/lrgsglib
      git submodule init
      git submodule update

2. **Create the conda environment with custom name:**

   If your project uses a different environment name:

   .. code-block:: bash

      conda env create -f lrgsglib/lrgsgenv.yml -n your_env_name
      conda activate your_env_name

3. **Point output paths at the outer project (optional):**

   By default the library derives its data/log/notebook directories from the
   ``lrgsglib`` subdirectory. To place them under the outer project root
   instead, set environment variables (or a ``.env`` file) before importing
   the library — from the outer project root:

   .. code-block:: bash

      export LRGSG_LLIB="$(pwd)"        # outer project root
      export LRGSG_DATA="$(pwd)/data"   # e.g. outer-project/data

   These are read by ``src/lrgsglib/config/lrgsg_env.py`` (environment first,
   otherwise derived from the package location). ``LRGSG_LOG`` and
   ``LRGSG_IPYNB`` are available too; all such directories are created lazily
   on first write.

4. **Install in editable mode:**

   From the outer project root:

   .. code-block:: bash

      pip install -e lrgsglib

Verifying Installation
----------------------

To verify the installation was successful:

.. code-block:: python

   import lrgsglib
   from lrgsglib import Lattice2D, ContactProcess
   import numpy as np

   # Create a small 2D lattice
   lattice = Lattice2D(side1=10, geo='sqr', pflip=0.2)
   print(f"Created lattice with {lattice.N} nodes")

   # Check C extensions are available
   cp = ContactProcess(lattice, gamma=0.5, runlang='C1c')
   print("C extensions loaded successfully!")

Running Tests
^^^^^^^^^^^^^

Run the test suite to ensure everything is working:

.. code-block:: bash

   # Full test suite
   pytest test/

   # Quick smoke test (small systems, minimal iterations)
   pytest test/ --quick

Building from Scratch
---------------------

To rebuild the compiled extensions after editing the C/C++ sources, re-run the
editable install — ``scikit-build-core`` recompiles the changed CMake/pybind11
targets:

.. code-block:: bash

   pip install -e .

For a submodule checkout, run the same command against the subdirectory from
the outer project root:

.. code-block:: bash

   pip install -e lrgsglib

Maintainers using the legacy ``make`` route can instead rebuild the standalone
C/C++ simulators directly (not required for normal use):

.. code-block:: bash

   make clean
   make all

Troubleshooting
---------------

**Common Issues:**

1. **Import errors after installation:**

   Ensure you've activated the correct conda environment and installed in editable mode.

2. **C extension compilation errors:**

   Check that you have a working C/C++ compiler:

   .. code-block:: bash

      gcc --version
      g++ --version

3. **Missing graph-tool:**

   graph-tool must be installed via conda, not pip:

   .. code-block:: bash

      conda install -c conda-forge graph-tool

4. **Path configuration issues in submodule mode:**

   Output paths are derived from the package location by default. To redirect
   them at the outer project, set the ``LRGSG_LLIB``/``LRGSG_DATA`` environment
   variables (or add an optional ``.env`` file) before importing the library,
   as shown in the Submodule Installation section above.

5. **Build system errors:**

   Inspect the ``pip install -e .`` output for the failing CMake/compiler
   step. Common issues include:

   - Missing build requirements (a C/C++ compiler or CMake)
   - Incorrect Python version (needs 3.12+)
   - Missing conda environment

Optional: GPU Support
---------------------

For GPU-accelerated operations using CuPy:

.. code-block:: bash

   conda install -c conda-forge cupy

Verify GPU support:

.. code-block:: python

   import cupy as cp
   print(f"CuPy version: {cp.__version__}")
   print(f"CUDA available: {cp.cuda.is_available()}")

Next Steps
----------

- :doc:`quickstart` - Quick tutorial to get started
- :doc:`user_guide/index` - Comprehensive user guide
- :doc:`api/index` - API reference documentation
