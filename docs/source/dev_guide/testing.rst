Testing
=======

This document describes test organization, how to run tests, write new tests,
and use benchmarks in lrgsglib.

Test Organization
-----------------

All tests live in the ``test/`` directory with specific naming conventions:

.. code-block:: text

   test/
   ├── .tmp/                           ← Temporary output files (gitignored)
   │   ├── *.pkl                       ← Benchmark results
   │   ├── *.png                       ← Generated plots
   │   └── *.log                       ← Log files
   │
   ├── conftest.py                     ← Shared fixtures + --quick/--show-plots
   ├── code_tests/                     ← Imports, instantiation, execution, portability
   ├── dynamics_tests/                 ← Ising, contact process, voter, Kuramoto, …
   ├── graph_tests/                    ← Graph construction, disorder, engine parity
   ├── integration_tests/              ← End-to-end pipelines
   ├── physical_tests/                 ← Physics / invariant checks
   ├── spectral_tests/                 ← Spectra, entropy, quantum propagator
   ├── utils_tests/                    ← Utility helpers
   ├── benchmarks/                     ← Performance benchmarks (bench_*.py)
   └── scripts/                        ← One-off utility scripts

Unit tests are ``test_*.py`` files inside the ``*_tests/`` domain
subdirectories.

File Naming Conventions
~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1

   * - Pattern
     - Purpose
     - Example
   * - ``test_*.py``
     - Unit tests (pytest)
     - ``test_contact_process.py``
   * - ``bench_*.py``
     - Performance benchmarks
     - ``bench_mcg_slaplspect_backends.py``
   * - ``demo_*.py``
     - Feature demonstrations
     - ``demo_edge_sign_clustering.py``
   * - ``diagnose_*.py``
     - Debugging specific issues
     - ``diagnose_cp_p0_anomalies.py``
   * - ``script_*.py``
     - One-off utilities
     - ``script_check_border_mode.py``

**Naming rules:**

- Be specific: ``bench_mcg_slaplspect_backends.py`` not ``benchmark.py``
- Include what's being tested: ``test_contact_process.py``
- Use descriptive names that explain the purpose

Running Tests
-------------

Quick Start
~~~~~~~~~~~

.. code-block:: bash

   # Run all tests
   cd lrgsglib
   pytest test/

   # Run with verbose output
   pytest test/ -v

   # Run a specific test file
   pytest test/dynamics_tests/test_contact_process.py

   # Run a specific test class (or ::Class::method for one function)
   pytest test/dynamics_tests/test_contact_process.py::TestContactProcess

   # Run tests matching a pattern
   pytest test/ -k "lattice"

Common pytest Options
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Show print statements
   pytest test/ -s

   # Stop on first failure
   pytest test/ -x

   # Run only previously failed tests
   pytest test/ --lf

   # Parallel execution
   pytest test/ -n auto

   # Coverage report
   pytest test/ --cov=src/lrgsglib --cov-report=html

   # Fast smoke run — small systems (side=4, steps=50); this is what CI runs
   pytest test/ --quick

   # Display generated figures interactively instead of only saving them
   pytest test/ --show-plots

.. note::

   ``--quick`` and ``--show-plots`` are custom options defined in
   ``test/conftest.py``; they drive the ``quick_mode``, ``system_size``,
   ``n_steps`` and ``save_figure`` fixtures. ``--quick`` is the flag the CI
   pipeline uses (see ``.github/workflows/quality.yml``).

Current Test Suite
------------------

Unit Tests
~~~~~~~~~~

.. code-block:: bash

   test/dynamics_tests/test_contact_process.py           # Contact process dynamics
   test/dynamics_tests/test_contact_process_limiting.py  # Edge cases
   test/spectral_tests/test_entropy_expm_multiply.py     # Matrix exponential
   test/spectral_tests/test_lattice2d_spectral.py        # Lattice spectral analysis
   test/dynamics_tests/test_lattice2d_animation.py       # Animation functionality
   test/integration_tests/test_mcg_output_format.py      # MCG output verification
   test/spectral_tests/test_quantum_propagator.py        # Quantum propagator
   test/dynamics_tests/test_voter_model.py               # Voter model dynamics
   test/graph_tests/test_refactoring.py                  # Refactoring checks

Benchmarks
~~~~~~~~~~

.. code-block:: bash

   # Contact process benchmarks
   test/benchmarks/bench_contact_process_quick.py     # Quick benchmark
   test/benchmarks/bench_contact_process_extended.py  # Extended benchmark
   test/benchmarks/bench_contact_process_c1b_vs_c1c.py  # Backend comparison

   # Spectral analysis benchmarks
   test/benchmarks/bench_mcg_slaplspect_backends.py   # Backend comparison
   test/benchmarks/bench_mcg_backends.py              # Alternative backends
   test/benchmarks/bench_entropy_methods.py           # Entropy computation methods

Writing Tests
-------------

Basic Test Structure
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   """Test module for contact process dynamics."""

   import pytest
   import numpy as np
   from lrgsglib.statsys.ContactProcess import ContactProcessEI
   from lrgsglib.graphs import Lattice2D


   @pytest.fixture
   def lattice():
       """Create a test lattice."""
       lat = Lattice2D(16, pflip=0.1, seed=42)
       lat.flip_random_fract_edges()
       return lat


   def test_contact_process_initialization(lattice):
       """Test ContactProcessEI initializes correctly."""
       cp = ContactProcessEI(sg=lattice, gamma=1.5, steps=100)
       assert cp.N == lattice.N
       assert cp.gamma == 1.5


   def test_contact_process_run(lattice):
       """Test ContactProcessEI completes a run."""
       cp = ContactProcessEI(sg=lattice, gamma=1.5, steps=100, runlang='py',
                             savedensity=True)
       cp.init_contact_dynamics()
       cp.run(verbose=False)
       assert len(cp.density) > 0


   @pytest.mark.parametrize("gamma", [0.5, 1.0, 2.0, 3.0])
   def test_contact_process_various_gamma(lattice, gamma):
       """Test ContactProcessEI with various gamma values."""
       cp = ContactProcessEI(sg=lattice, gamma=gamma, steps=50, runlang='py')
       cp.init_contact_dynamics()
       cp.run(verbose=False)
       # Should complete without error

Test Fixtures
~~~~~~~~~~~~~

Use fixtures for common setup:

.. code-block:: python

   import pytest
   from pathlib import Path

   # Temporary directory for test outputs
   TEST_TMP_DIR = Path(__file__).parent / ".tmp"


   @pytest.fixture(scope="module")
   def test_tmp_dir():
       """Ensure .tmp directory exists."""
       TEST_TMP_DIR.mkdir(exist_ok=True)
       return TEST_TMP_DIR


   @pytest.fixture
   def small_lattice():
       """Small lattice for fast tests."""
       lat = Lattice2D(8, pflip=0.2, seed=123)
       lat.flip_random_fract_edges()
       return lat


   @pytest.fixture
   def medium_lattice():
       """Medium lattice for integration tests."""
       lat = Lattice2D(32, pflip=0.3, seed=456)
       lat.flip_random_fract_edges()
       return lat

Parametrized Tests
~~~~~~~~~~~~~~~~~~

Test multiple inputs efficiently:

.. code-block:: python

   @pytest.mark.parametrize("side,pflip", [
       (8, 0.0),
       (8, 0.1),
       (16, 0.2),
       (32, 0.3),
   ])
   def test_lattice_creation(side, pflip):
       """Test lattice creation with various parameters."""
       lat = Lattice2D(side, pflip=pflip, seed=42)
       assert lat.N == side * side


   @pytest.mark.parametrize("backend", ["numpy", "scipy"])
   def test_spectrum_backends(small_lattice, backend):
       """Test spectral computation with different backends."""
       small_lattice.compute_laplacian_spectrum(backend=backend)
       assert small_lattice.eigv is not None

Markers and Skipping
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import pytest

   @pytest.mark.slow
   def test_large_computation():
       """Test that takes a long time."""
       pass

   @pytest.mark.skipif(
       not has_cupy(),
       reason="CuPy not available"
   )
   def test_gpu_backend():
       """Test GPU computation."""
       pass

   @pytest.mark.xfail(reason="Known issue #123")
   def test_known_failure():
       """Test with known failure."""
       pass

Running marked tests:

.. code-block:: bash

   # Skip slow tests
   pytest test/ -m "not slow"

   # Run only slow tests
   pytest test/ -m slow

Writing Benchmarks
------------------

Benchmark Structure
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   #!/usr/bin/env python
   """
   Benchmark MCG_SlaplSpect.py with different backends.

   Compares scipy vs cupy performance across graph sizes.
   Output files stored in test/.tmp/
   """

   from pathlib import Path
   import pickle
   import time
   import numpy as np

   # Output directory
   TEST_TMP_DIR = Path(__file__).parent / ".tmp"
   TEST_TMP_DIR.mkdir(exist_ok=True)


   def benchmark_scipy(sizes):
       """Benchmark scipy backend."""
       results = []
       for N in sizes:
           start = time.time()
           # ... computation ...
           elapsed = time.time() - start
           results.append({"N": N, "time": elapsed, "backend": "scipy"})
       return results


   def benchmark_cupy(sizes):
       """Benchmark cupy backend."""
       try:
           import cupy
       except ImportError:
           return []
       results = []
       for N in sizes:
           start = time.time()
           # ... computation ...
           elapsed = time.time() - start
           results.append({"N": N, "time": elapsed, "backend": "cupy"})
       return results


   def main():
       sizes = [100, 500, 1000, 2000, 5000]

       results = []
       results.extend(benchmark_scipy(sizes))
       results.extend(benchmark_cupy(sizes))

       # Save results
       output_file = TEST_TMP_DIR / "bench_backends_results.pkl"
       with open(output_file, "wb") as f:
           pickle.dump(results, f)

       print(f"Results saved to {output_file}")


   if __name__ == "__main__":
       main()

Output File Management
----------------------

**Important:** Keep the test folder clean!

1. **Temporary files** go to ``test/.tmp/``:

   .. code-block:: python

      TEST_TMP_DIR = Path(__file__).parent / ".tmp"
      TEST_TMP_DIR.mkdir(exist_ok=True)
      output_file = TEST_TMP_DIR / "results.pkl"

2. **Documentation** goes to ``.agents/devd/<project>/``:

   - Never put markdown files in the test folder
   - Agent reports belong in ``.agents/``

3. **Persistent outputs** (if necessary) go to ``test/output/`` or ``test/figures/``:

   - Use sparingly
   - Only for outputs that need to be committed

Testing Best Practices
----------------------

**Reproducibility:**

- Always use fixed seeds for random operations
- Document expected outputs in test names

**Speed:**

- Keep unit tests fast (< 1 second each)
- Mark slow tests with ``@pytest.mark.slow``
- Use small test data

**Isolation:**

- Tests should not depend on each other
- Clean up temporary files
- Use fixtures for setup/teardown

**Coverage:**

- Test edge cases and error conditions
- Test with various parameter combinations
- Test both Python and C backends

Doctests
--------

lrgsglib uses doctests for inline documentation testing:

.. code-block:: bash

   # Run doctests in documentation
   make -C docs doctest

   # Run doctests in source files
   pytest --doctest-modules src/lrgsglib/

Writing doctests:

.. code-block:: python

   def compute_mean(values):
       """
       Compute arithmetic mean.

       Parameters
       ----------
       values : array-like
           Input values.

       Returns
       -------
       float
           The mean.

       Examples
       --------
       >>> compute_mean([1, 2, 3, 4, 5])
       3.0
       >>> compute_mean([10])
       10.0
       """
       return sum(values) / len(values)

Continuous Integration
----------------------

CI is driven by three GitHub Actions workflows under
``.github/workflows/``:

**Code Quality** (``quality.yml``) — runs on pushes to ``main``,
``dev-notebooks`` and ``feature/**``, and on pull requests to ``main`` and
``dev-notebooks``:

1. **Lint** — ``black --check``, ``isort --check`` and ``flake8`` over
   ``src/lrgsglib/``.
2. **Docs** — build the Sphinx site with ``make html`` under ``-W``
   (warnings are errors).
3. **Test** — ``pip install -e ".[dev]"`` then ``pytest test/ --quick`` on
   Ubuntu and macOS.

**Build & Test Wheels** (``build-wheels.yml``) — runs on pushes to
``main``/``dev``/``feature/**``, on pull requests to ``main``/``dev``, and on
``v*`` tags. It builds wheels with ``cibuildwheel`` and an sdist, installs the
built wheel and runs ``pytest test/ --quick`` against it, then publishes to
PyPI on tag pushes.

**Deploy Documentation** (``docs-deploy.yml``) — runs on pushes to ``main``
and on ``v*`` tags. It builds the docs with micromamba and deploys the
multi-version site to the ``gh-pages`` branch.

.. note::

   Tests install via ``pip`` (``.[dev]``), not conda; C extensions are
   compiled implicitly by the editable / ``cibuildwheel`` install. Conda
   (micromamba) is used only by the documentation-deploy workflow.

See Also
--------

- :doc:`contributing` - Development workflow
- :doc:`style_guide` - Code conventions
- :doc:`build_system` - Building for tests
