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
- pybind11 >= 2.10.0
- joblib >= 1.0.0
- numba >= 0.58.0

Optional dependencies:

- cupy >= 9.0.0 (for GPU acceleration)
- graph-tool (from conda, for high-performance graph operations)

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

3. **Build the project:**

   This step compiles C/C++ extensions and configures environment variables:

   .. code-block:: bash

      make all

   The build process:

   - Compiles C/C++ simulators in ``src/lrgsglib/Ccore/``
   - Builds pybind11 Python extensions
   - Generates a ``.env`` file with configured paths
   - Sets up data, log, and binary directories

4. **Install in editable mode:**

   .. code-block:: bash

      pip install -e .

   For development with additional tools:

   .. code-block:: bash

      pip install -e ".[dev,jupyter,docs]"

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

3. **Build with custom paths:**

   Configure paths relative to the outer project root instead of the ``lrgsglib`` subdirectory:

   .. code-block:: bash

      cd lrgsglib
      make all LRGSG_LLIB=$(pwd)/.. CONDA_ENV_NAME=your_env_name

   **What this does:**

   - Sets ``LRGSG_LLIB`` to the outer project root
   - Configures data folder as ``outer-project/data`` (not ``lrgsglib/data``)
   - Configures notebooks as ``outer-project/ipynb`` (not ``lrgsglib/ipynb``)
   - Configures logs as ``outer-project/.log`` (not ``lrgsglib/.log``)
   - Keeps library source code paths relative to ``lrgsglib/``
   - Uses your custom conda environment name

   **Example for lrgsglib-ipynb:**

   .. code-block:: bash

      cd lrgsglib
      make all LRGSG_LLIB=$(pwd)/.. CONDA_ENV_NAME=lrgsgnb

4. **Install in editable mode:**

   .. code-block:: bash

      cd ..  # Return to outer project root
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

   # Quick smoke test
   python test/quick_test.py

   # Extended validation
   python test/extended_test.py

Building from Scratch
---------------------

If you need to rebuild after making changes:

.. code-block:: bash

   make clean
   make all

For submodule builds, include the same parameters:

.. code-block:: bash

   make clean
   make all LRGSG_LLIB=$(pwd)/.. CONDA_ENV_NAME=your_env_name

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

   Verify the ``.env`` file in ``lrgsglib/`` contains correct paths pointing to your outer project directories.

5. **Build system errors:**

   Check the Makefile output for specific error messages. Common issues include:

   - Missing dependencies (numpy, pybind11)
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
