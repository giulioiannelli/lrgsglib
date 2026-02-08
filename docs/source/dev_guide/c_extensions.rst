C/C++ Extensions
================

This document describes the C/C++ performance layer in lrgsglib, including
the simulator binaries, random number generation, and Python integration.

Overview
--------

Performance-critical simulations are implemented in C for ~100x speedup over
pure Python. C code is co-located with Python classes in a folder-per-class
layout under ``src/lrgsglib/statsys/``. Each model has its own ``ccore/``
subdirectory, and shared infrastructure lives in ``statsys/_ccore/``.

Directory Structure
-------------------

.. code-block:: text

   statsys/
   ├── _ccore/                      ← Shared C infrastructure
   │   ├── SFMT/                    ← SIMD-oriented Fast Mersenne Twister
   │   │   ├── SFMT.c/.h            ← Core RNG implementation
   │   │   └── params/              ← Period parameters (19937, etc.)
   │   ├── LRGSG_utils.c/.h         ← Common utility functions
   │   ├── LRGSG_customs.h          ← Custom type definitions
   │   ├── LRGSG_bindynsys.c/.h     ← Binary dynamics utilities
   │   ├── LRGSG_contdynsys.c/.h    ← Continuous dynamics utilities
   │   ├── LRGSG_vecdynsys.c/.h     ← Vector dynamics utilities
   │   └── sfmtrng.c/.h             ← SFMT RNG wrapper
   │
   ├── IsingDynamics/               ← Ising model
   │   ├── IsingDynamics.py
   │   └── ccore/
   │       ├── bin/                  ← Compiled binaries (IsingSimulator*)
   │       ├── LRGSG_rbim.c/.h
   │       ├── base/                 ← Core Ising functions (pybind11)
   │       └── storer/               ← Data storage (pybind11)
   │
   ├── ContactProcess/ccore/        ← Contact process simulators
   ├── VoterModel/ccore/            ← Voter model simulators
   ├── KuramotoModel/ccore/         ← Kuramoto oscillators
   ├── PottsModel/ccore/            ← Potts model
   ├── XYModel/ccore/               ← XY (planar rotator) model
   ├── HeisenbergModel/ccore/       ← Heisenberg model
   ├── MultiSpeciesModel/ccore/     ← Multi-species model
   └── SignedRW/ccore/              ← Signed random walks

Simulator Variants
------------------

Ising Model Simulators
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1

   * - Binary
     - Algorithm
     - Features
     - ``runlang``
   * - ``IsingSimulator0``
     - Basic Metropolis
     - Minimal, fast
     - ``C0``
   * - ``IsingSimulator1b``
     - Optimized Metropolis
     - Logging, magnetization tracking
     - ``C1b``
   * - ``IsingSimulator3b``
     - Simulated Annealing
     - Temperature schedule
     - ``C3b``
   * - ``IsingSimulator4b``
     - Parallel Tempering
     - Replica exchange
     - ``C4b``

**Usage from Python:**

.. code-block:: python

   from lrgsglib.statsys.IsingDynamics import IsingDynamics
   from lrgsglib.graphs.Lattice2D import Lattice2D

   lattice = Lattice2D(side=64, pflip=0.3)
   lattice.flip_random_fract_edges()

   # Use C backend with runlang parameter
   ising = IsingDynamics(sg=lattice, T=2.0, steps=10000, runlang='C1b')
   ising.init_ising_dynamics()
   ising.run(verbose=False)

Contact Process Simulators
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1

   * - Binary
     - Description
     - ``runlang``
   * - ``ContactSimulator0``
     - Basic implementation
     - ``C0``
   * - ``ContactSimulator1c``
     - Optimized with survival tracking
     - ``C1c``
   * - ``ContactSimulator1d``
     - 1D lattice optimized
     - ``C1d``

**Usage from Python:**

.. code-block:: python

   from lrgsglib.statsys.ContactProcess import ContactProcessEI

   cp = ContactProcessEI(sg=lattice, gamma=1.5, steps=5000, runlang='C1c')
   cp.init_contact_dynamics()
   cp.run(verbose=False)

Voter Model Simulators
~~~~~~~~~~~~~~~~~~~~~~

- ``VoterSimulator0`` - Basic voter dynamics
- ``VoterSimulator1`` - With magnetization logging

Random Number Generation
------------------------

All C simulators use the SIMD-oriented Fast Mersenne Twister (SFMT) for
high-quality, fast random number generation.

SFMT Features
~~~~~~~~~~~~~

- Period: 2^19937 - 1 (default)
- SIMD vectorized for modern CPUs
- Uniform distribution in [0, 1)
- Jump-ahead capability for parallel streams

**Initialization:**

.. code-block:: c

   #include "sfmtrng.h"

   // Initialize with seed
   sfmt_t sfmt;
   sfmt_init_gen_rand(&sfmt, seed);

   // Generate random numbers
   double r = sfmt_genrand_real2(&sfmt);  // [0, 1)
   uint32_t i = sfmt_genrand_uint32(&sfmt);

Common C Utilities
------------------

LRGSG_utils.c
~~~~~~~~~~~~~

General-purpose utilities used across all simulators:

.. code-block:: c

   // File I/O
   void read_binary_int8(const char* path, int8_t* data, int N);
   void write_binary_int8(const char* path, int8_t* data, int N);

   // Array operations
   void init_array_int8(int8_t* arr, int N, int8_t val);
   double compute_mean_int8(int8_t* arr, int N);

LRGSG_bindynsys.c
~~~~~~~~~~~~~~~~~

Binary dynamics system utilities:

.. code-block:: c

   // Edge list management
   typedef struct {
       int* src;
       int* dst;
       int8_t* sign;
       int E;
   } EdgeList;

   void load_edge_list(const char* path, EdgeList* el);
   void free_edge_list(EdgeList* el);

Ising-Specific Code
~~~~~~~~~~~~~~~~~~~

Located in ``statsys/RBIsingM/``:

.. code-block:: c

   // LRGSG_rbim.h - Random-bond Ising model
   double compute_energy(int8_t* s, EdgeList* el, int N);
   double compute_magnetization(int8_t* s, int N);
   void metropolis_step(int8_t* s, EdgeList* el, double beta, sfmt_t* rng);

   // LRGSG_sa.h - Simulated annealing
   void simulated_annealing(int8_t* s, EdgeList* el, double T_start,
                            double T_end, int steps, sfmt_t* rng);

   // LRGSG_pt.h - Parallel tempering
   void parallel_tempering(int8_t** replicas, EdgeList* el,
                           double* temps, int n_replicas, int steps,
                           sfmt_t* rng);

Python-C Communication
----------------------

Data Exchange Protocol
~~~~~~~~~~~~~~~~~~~~~~

Python and C communicate through binary files and stdout:

**Input (Python → C):**

1. **Edge list** - Binary file with graph connectivity and signs
2. **Adjacency** - Binary adjacency list for efficient neighbor lookup
3. **Initial state** - Binary file with initial spin configuration

**Output (C → Python):**

1. **Final state** - Via stdout (binary int8 array)
2. **Observables** - Binary files (magnetization, energy time series)

File Formats
~~~~~~~~~~~~

**Edge list format** (``.edgl.bin``):

.. code-block:: text

   [int32: E (number of edges)]
   [int32 * E: source nodes]
   [int32 * E: destination nodes]
   [int8 * E: edge signs (+1 or -1)]

**Spin configuration** (``.s.bin``):

.. code-block:: text

   [int8 * N: spin values (+1 or -1)]

**Magnetization** (``.m.bin``):

.. code-block:: text

   [float64 * T: magnetization at each time step]

C Program Interface
~~~~~~~~~~~~~~~~~~~

All simulators follow this command-line interface:

.. code-block:: bash

   ./IsingSimulator1b N pflip steps datadir syshape run_id out_id

Arguments:

- ``N`` - Number of nodes
- ``pflip`` - Fraction of negative edges
- ``steps`` - Number of Monte Carlo sweeps
- ``datadir`` - Directory containing input files
- ``syshape`` - System shape descriptor (e.g., "L=64")
- ``run_id`` - Suffix for input files
- ``out_id`` - Suffix for output files

Building C Extensions
---------------------

The Makefile handles compilation:

.. code-block:: bash

   # Build all C programs
   make c-make

   # Build specific simulator
   make src/lrgsglib/statsys/IsingDynamics/ccore/bin/IsingSimulator1b

   # Clean and rebuild
   make clean && make c-make

Compilation flags (from ``build/cconfig.mk``):

.. code-block:: make

   GCC = gcc
   CFLAGS = -O3 -march=native -Wall
   LMFLAG = -lm
   ALLFLAGS = $(CFLAGS) -DSFMT_MEXP=19937

Adding New Simulators
---------------------

To add a new simulator:

1. **Create C source file:**

   .. code-block:: c

      // statsys/MyModel/ccore/MySimulator.c
      #include "LRGSG_utils.h"    /* resolved via -I flags */
      #include "sfmtrng.h"
      #include "LRGSG_mymodel.h"

      int main(int argc, char* argv[]) {
          // Parse arguments
          int N = atoi(argv[1]);
          // ... simulation logic ...
          // Write final state to stdout
          fwrite(state, sizeof(int8_t), N, stdout);
          return 0;
      }

2. **Add Makefile rule:**

   .. code-block:: make

      $(LRGSG_CCORE_BIN)/MySimulator%: $(LRGSG_STATSYS_MYMODEL)/MySimulator%.c \
                                        $(PATH_SRCC_FILES) \
                                        $(PATH_SFMT_FILES)
          $(GCC) $(ALLFLAGS) -o $@ $^ $(LMFLAG)

3. **Create Python wrapper:**

   .. code-block:: python

      class MyModel(BinDynSys):
          def build_cprogram_command(self):
              self.CbaseName = f"MySimulator{self.runlang[-1]}"
              self.cprogram = [LRGSG_CCORE_BIN / self.CbaseName] + args

          def run_cprogram(self, verbose=False):
              result = subprocess.run(self.cprogram, ...)
              self.s = np.frombuffer(result.stdout, dtype=np.int8)

pybind11 Bindings
-----------------

Some C++ code uses pybind11 for direct Python bindings (in ``bindings/``):

.. code-block:: cpp

   // bindings/random_walk.cpp
   #include <pybind11/pybind11.h>
   #include <pybind11/numpy.h>

   namespace py = pybind11;

   py::array_t<double> random_walk(int steps, double p) {
       // ... implementation ...
   }

   PYBIND11_MODULE(random_walk, m) {
       m.def("random_walk", &random_walk, "Perform random walk");
   }

**Building pybind11 extensions:**

.. code-block:: bash

   # Compiled automatically with pip install -e .
   pip install -e .

Performance Considerations
--------------------------

**Memory:**

- C simulators use compact int8 for spins (1 byte vs 8 bytes for Python int)
- Edge lists stored as contiguous arrays for cache efficiency

**Computation:**

- SIMD-vectorized RNG reduces random number overhead
- Locality-aware algorithms minimize cache misses
- No Python GIL contention during simulation

**I/O:**

- Binary files are ~10x faster than text for large arrays
- stdout used for small outputs to avoid file system overhead

Debugging C Code
----------------

**Enable debug symbols:**

.. code-block:: bash

   CFLAGS="-g -O0" make c-make

**Run with valgrind:**

.. code-block:: bash

   valgrind --leak-check=full ./src/lrgsglib/statsys/IsingDynamics/ccore/bin/IsingSimulator1b ...

**Add debug prints:**

.. code-block:: c

   #ifdef DEBUG
   fprintf(stderr, "Step %d: E=%.3f M=%.3f\n", step, energy, magn);
   #endif

See Also
--------

- :doc:`architecture` - System overview
- :doc:`build_system` - Build process details
- :doc:`testing` - Testing C code
