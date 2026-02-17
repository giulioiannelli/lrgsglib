Statistical Physics Simulations
================================

This guide covers the statistical physics simulation capabilities of lrgsglib,
including Ising model dynamics, contact processes, and voter models on signed
graphs.

.. contents:: Contents
   :local:
   :depth: 2

Overview
--------

lrgsglib provides three core dynamical systems:

1. **Ising Model** (``IsingDynamics``): Spin dynamics with Metropolis, Simulated
   Annealing (SA), and Parallel Tempering (PT)
2. **Contact Process** (``ContactProcessSIR``, ``ContactProcessEI``): Infection
   spreading and excitation-inhibition dynamics
3. **Voter Model** (``VoterModel``): Opinion dynamics with neighbor copying

All dynamics classes share a common pattern:

.. code-block:: python

   # 1. Create a signed graph
   graph = Lattice2D(side1=32, pflip=0.3, seed=42)
   graph.flip_random_fract_edges()

   # 2. Initialize dynamics (use sg= parameter!)
   dynamics = IsingDynamics(sg=graph, T=2.0, steps=1000, runlang='C1b')

   # 3. Call the init method
   dynamics.init_ising_dynamics()

   # 4. Run the simulation
   dynamics.run(verbose=False)

   # 5. Access results
   print(dynamics.magn[-1])  # Final magnetization

.. important::

   Always use ``sg=`` (not ``graph=``) when passing the signed graph to dynamics classes.
   This is a common source of errors!

Ising Model Dynamics
--------------------

The Ising model simulates spin systems where each node has a spin :math:`s_i \in \{-1, +1\}`.
On signed graphs, positive edges favor alignment while negative edges favor anti-alignment.

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   from lrgsglib.graphs import Lattice2D
   from lrgsglib.statsys.IsingDynamics import IsingDynamics

   # Create a frustrated lattice
   lattice = Lattice2D(side1=32, geo='sqr', pflip=0.3, seed=42)
   lattice.flip_random_fract_edges()

   # Initialize Ising dynamics
   ising = IsingDynamics(
       sg=lattice,          # SignedGraph instance (NOT graph=!)
       T=2.0,               # Temperature (NOT temperature=!)
       steps=10000,         # Number of Monte Carlo sweeps
       runlang='C1b',       # Use fast C backend
   )

   # Initialize and run
   ising.init_ising_dynamics()
   ising.run(verbose=False)

   # Access results
   print(f"Final magnetization: {ising.magn[-1]:.4f}")
   print(f"Final energy: {ising.ene[-1]:.4f}")
   print(f"Number of data points: {len(ising.magn)}")

Backend Selection
~~~~~~~~~~~~~~~~~

The ``runlang`` parameter selects the execution backend and algorithm.
The canonical format is ``<backend>_<algorithm>[_<output>]``:

.. list-table:: Backend Identifiers
   :header-rows: 1
   :widths: 10 15 75

   * - Backend
     - Prefix
     - Description
   * - Python
     - ``py_``
     - Pure-Python Metropolis / SA / PT. Always available, slowest.
   * - C subprocess
     - ``c_``
     - Compiled C binaries via file I/O. ~100x faster. **NX graphs only.**
   * - Pybind11
     - ``pb_``
     - C kernels called in-process via numpy arrays. No file I/O. Works with any graph.
   * - CuPy (GPU)
     - ``cu_``
     - CUDA kernels on GPU via CuPy. Fastest for large systems. Works with any graph.

.. list-table:: Algorithm Identifiers
   :header-rows: 1
   :widths: 10 15 75

   * - Algorithm
     - Suffix
     - Description
   * - Metropolis
     - ``_met``
     - Glauber-Metropolis at fixed temperature
   * - Simulated Annealing
     - ``_sa``
     - Temperature cooling schedule
   * - Parallel Tempering
     - ``_pt``
     - Replica exchange at multiple temperatures
   * - Wolff Cluster
     - ``wolff``
     - Single-cluster update (no prefix needed)
   * - Swendsen-Wang
     - ``sw``
     - Multi-cluster update (no prefix needed)

**Examples:** ``pb_met`` (pybind11 Metropolis), ``cu_sa`` (GPU simulated annealing),
``c_met_em`` (C subprocess Metropolis with E/M output), ``wolff`` (Wolff cluster).

**Convenience aliases:** ``"py"`` = ``py_met``, ``"pb"`` = ``pb_met``,
``"cu"`` / ``"cuda"`` / ``"gpu"`` = ``cu_met``.

**Legacy codes** (``"C1b"``, ``"C3b"``, ``"C4b"``) are mapped automatically
for backward compatibility.

.. list-table:: Backend Compatibility Matrix
   :header-rows: 1
   :widths: 20 15 15 15 15

   * - Graph Type
     - ``py_*``
     - ``c_*``
     - ``pb_*``
     - ``cu_*``
   * - NX (``Lattice2D``, etc.)
     - Yes
     - Yes
     - Yes
     - Yes
   * - GT (``Lattice2DGT``, etc.)
     - Yes
     - No (auto-fallback)
     - Yes
     - Yes

.. note::

   When a C subprocess backend is requested for a graph-tool graph, the system
   automatically falls back to the Python backend with a ``RuntimeWarning``.

Simulated Annealing
~~~~~~~~~~~~~~~~~~~

Find ground states via temperature cooling:

.. code-block:: python

   # Simulated Annealing configuration
   ising_sa = IsingDynamics(
       sg=lattice,
       sa_enabled=True,
       T_init=10.0,              # Starting temperature
       T_final=0.01,             # Final temperature
       cooling_schedule='exponential',  # 'linear', 'exponential', 'logarithmic'
       cooling_rate=0.95,        # Factor for exponential cooling
       steps_per_T=100,          # MC sweeps per temperature
       n_temperatures=100,       # Number of temperature steps
       runlang='C3b',
   )

   ising_sa.init_ising_dynamics()
   ising_sa.run(verbose=False)

   # Access cooling trajectory
   print(f"Final energy: {ising_sa.ene[-1]:.4f}")

Parallel Tempering
~~~~~~~~~~~~~~~~~~

Enhanced sampling with replica exchange:

.. code-block:: python

   # Parallel Tempering configuration
   ising_pt = IsingDynamics(
       sg=lattice,
       pt_enabled=True,
       n_replicas=8,             # Number of temperature replicas
       T_min=0.5,                # Lowest temperature
       T_max=5.0,                # Highest temperature
       T_ladder_type='geometric',  # 'geometric', 'linear', 'custom'
       steps_per_exchange=10,    # MC sweeps between exchange attempts
       n_exchanges=1000,         # Total number of exchange cycles
       runlang='C4b',
   )

   ising_pt.init_ising_dynamics()
   ising_pt.run(verbose=False)

Pybind11 Backend (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The pybind11 backend calls C kernels in-process via numpy arrays.
No file I/O is needed, and it works with **any graph type** (NX or GT):

.. code-block:: python

   ising_pb = IsingDynamics(
       sg=lattice,
       T=2.0,
       steps=1000,
       runlang='pb_met',  # Pybind11 Metropolis
       seed=42,
   )
   ising_pb.init_ising_dynamics()
   ising_pb.run(verbose=False)

Works with graph-tool graphs too:

.. code-block:: python

   from lrgsglib.graphs.gt.Lattice2DGT import Lattice2DGT

   gt_lattice = Lattice2DGT(side1=32, geo='sqr', pflip=0.2, seed=42)
   gt_lattice.flip_random_fract_edges()

   ising_gt = IsingDynamics(sg=gt_lattice, T=2.0, steps=1000, runlang='pb_met')
   ising_gt.init_ising_dynamics()
   ising_gt.run(verbose=False)

CuPy GPU Backend
~~~~~~~~~~~~~~~~

For large systems, the GPU backend provides the fastest execution.
Requires CuPy and a CUDA-capable GPU:

.. code-block:: python

   ising_gpu = IsingDynamics(
       sg=lattice,
       T=2.0,
       steps=5000,
       runlang='cu_met',  # GPU Metropolis (aliases: 'cuda', 'gpu')
       seed=42,
   )
   ising_gpu.init_ising_dynamics()
   ising_gpu.run(verbose=False)

The GPU backend uses a graph-coloring approach: nodes are partitioned into
independent sets (colors) so all same-color nodes can be updated simultaneously.
For bipartite graphs (e.g. square lattices), this is the classic checkerboard
decomposition with 2 colors. Works on arbitrary graph topologies.

GPU Simulated Annealing and Parallel Tempering are also available:

.. code-block:: python

   # GPU SA
   ising_sa = IsingDynamics(sg=lattice, runlang='cu_sa', sa_enabled=True,
       T_init=5.0, T_final=0.01, n_temperatures=100)

   # GPU PT
   ising_pt = IsingDynamics(sg=lattice, runlang='cu_pt', pt_enabled=True,
       n_replicas=8, T_min=0.5, T_max=5.0, n_exchanges=500)

Python Backend
~~~~~~~~~~~~~~

For debugging or when compiled backends aren't available:

.. code-block:: python

   ising_py = IsingDynamics(
       sg=lattice,
       T=2.0,
       steps=1000,
       runlang='py',  # Pure Python (slower but always available)
   )

   ising_py.init_ising_dynamics()
   ising_py.run(tqdm_on=True)  # Show progress bar

Cluster Algorithms (Wolff & Swendsen-Wang)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cluster algorithms avoid critical slowing down near phase transitions by
flipping entire clusters of aligned spins rather than single sites:

.. code-block:: python

   # Wolff algorithm — single-cluster update
   ising_wolff = IsingDynamics(
       sg=lattice, T=2.27, steps=5000,
       runlang='wolff',
   )
   ising_wolff.init_ising_dynamics()
   ising_wolff.run(verbose=False)

   # Swendsen-Wang — multi-cluster update (all clusters flipped per sweep)
   ising_sw = IsingDynamics(
       sg=lattice, T=2.27, steps=5000,
       runlang='sw',
   )
   ising_sw.init_ising_dynamics()
   ising_sw.run(verbose=False)

Cluster algorithms are pure-Python implementations that work with any graph
type (NX or GT). They are especially effective near the critical temperature
where single-spin Metropolis suffers from critical slowing down.

.. note::

   Cluster algorithms treat edge signs as bond weights. On frustrated
   (signed) graphs, they reduce to the standard Wolff/SW algorithm applied
   to the signed coupling :math:`J_{ij} \cdot s_i \cdot s_j`.

.. _topological-algorithms:

Topological Algorithms
~~~~~~~~~~~~~~~~~~~~~~

Topological algorithms exploit the **Laplacian spectral structure** of the graph
to guide Ising ground-state search. They are especially effective on frustrated
lattices where standard Metropolis gets trapped in metastable states.

**Topological Metropolis** (``topo_met``):
Monte Carlo in spectral coefficient space. Instead of flipping individual spins,
the algorithm perturbs coefficients in a decomposition over Laplacian
eigenvectors and binarizes the resulting continuous field:

.. code-block:: python

   # Topological Metropolis on a frustrated lattice
   ising_topo = IsingDynamics(
       sg=lattice,
       T=0.5,
       steps=10000,
       runlang='py_topo_met',   # or 'pb_topo_met' for pybind11
       topo_n_modes=40,         # Number of eigenvectors to use
       topo_polish=True,        # T=0 greedy polish per chunk
       topo_polish_sweeps=50,   # Max sweeps for polish
       topo_chunk_size=5000,    # Adaptive sigma adjustment interval
       topo_sigma_init=0.15,    # Initial proposal width
       seed=42,
   )
   ising_topo.init_ising_dynamics()
   ising_topo.run(verbose=False)

   # Results
   print(f"Best energy found: {ising_topo.topo_met_best_energy:.2f}")
   print(f"Final coefficients: {ising_topo.topo_met_coeffs}")
   best_config = ising_topo.topo_met_best_spins

Key features:

- **Mode-energy-weighted selection**: modes with lower RBIM energy are
  perturbed more frequently (:math:`P(m) \propto |E_{\text{RBIM}}[m]|`)
- **Adaptive sigma**: per-mode proposal widths are adjusted via Robbins-Monro
  to target ~30% acceptance rate
- **Greedy polish**: optional T=0 Metropolis sweeps after each chunk

**Topological Field Cooling Annealing** (``topo_fca``):
Standard simulated annealing with a static external field constructed from
softmax-weighted eigenvectors:

.. code-block:: python

   # TFCA on a frustrated lattice
   ising_tfca = IsingDynamics(
       sg=lattice,
       T=1.0,
       steps=100,
       runlang='py_topo_fca',   # or 'pb_topo_fca' for pybind11
       sa_enabled=True,
       T_init=5.0,
       T_final=0.01,
       n_temperatures=100,
       topo_n_modes=40,
       topo_tau=1.0,             # Softmax temperature
       topo_field_strength=1.0,  # Field magnitude scale
       seed=42,
   )
   ising_tfca.init_ising_dynamics()
   ising_tfca.run(verbose=False)

   # The spectral field is stored for analysis
   spectral_field = ising_tfca.topo_fca_field

The field is computed as:

.. math::

   h_i = \text{strength} \cdot \sum_k c_k \, v_k(i), \quad
   c_k = \text{softmax}(-E_k / \tau)

where :math:`E_k` is the RBIM energy of the *k*-th binarized eigenvector and
:math:`v_k` is the continuous eigenvector. Lower tau concentrates weight on
the best modes.

Both algorithms are available with Python (``py_``) and pybind11 (``pb_``)
backends. Aliases ``topo_met`` and ``topo_fca`` default to the Python backend.

Graph Engine Selection
~~~~~~~~~~~~~~~~~~~~~~

The ``--graph_engine`` (``-ge``) CLI flag selects between NetworkX (``nx``)
and graph-tool (``gt``) backends for standalone programs:

.. code-block:: bash

   # Default: NetworkX
   python src/L2D_IsingDynamics.py 16 0.0 -T 2.0 -rl py -na 1

   # graph-tool backend
   python src/L2D_IsingDynamics.py 16 0.0 -T 2.0 -rl py -na 1 -ge gt

When using the GT engine:

- **C subprocess backends** (``c_*``) are **not available** — the system
  automatically falls back to the Python backend with a warning.
- **Pybind11** (``pb_*``) and **CuPy** (``cu_*``) backends work with GT graphs.
- **Cluster algorithms** (``wolff``, ``sw``) work with GT graphs.
- GT graphs use a **neighbor cache** for dynamics (lazily built on first
  access), giving performance within ~1.3x of NX.

Complete Runlang Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Code
     - Backend
     - Algorithm
     - Notes
   * - ``py`` / ``py_met``
     - Python
     - Metropolis
     - Always available, slowest
   * - ``py_sa``
     - Python
     - Simulated Annealing
     - Always available
   * - ``py_pt``
     - Python
     - Parallel Tempering
     - Always available
   * - ``c_met`` / ``C1b``
     - C subprocess
     - Metropolis
     - NX only, ~100x faster
   * - ``c_sa`` / ``C3b``
     - C subprocess
     - Simulated Annealing
     - NX only
   * - ``c_pt`` / ``C4b``
     - C subprocess
     - Parallel Tempering
     - NX only
   * - ``pb_met`` / ``pb``
     - Pybind11
     - Metropolis
     - In-process C, any graph
   * - ``pb_sa``
     - Pybind11
     - Simulated Annealing
     - In-process C, any graph
   * - ``pb_pt``
     - Pybind11
     - Parallel Tempering
     - In-process C, any graph
   * - ``cu_met`` / ``cu`` / ``gpu``
     - CuPy (GPU)
     - Metropolis
     - Requires CUDA, any graph
   * - ``cu_sa``
     - CuPy (GPU)
     - Simulated Annealing
     - Requires CUDA, any graph
   * - ``cu_pt``
     - CuPy (GPU)
     - Parallel Tempering
     - Requires CUDA, any graph
   * - ``wolff``
     - Python
     - Wolff cluster
     - Any graph, near-Tc efficient
   * - ``sw``
     - Python
     - Swendsen-Wang cluster
     - Any graph, near-Tc efficient
   * - ``topo_met`` / ``py_topo_met``
     - Python
     - Topological Metropolis
     - Spectral coefficient space MC
   * - ``pb_topo_met``
     - Pybind11
     - Topological Metropolis
     - In-process C, any graph
   * - ``topo_fca`` / ``py_topo_fca``
     - Python
     - Topological FCA
     - SA with spectral field
   * - ``pb_topo_fca``
     - Pybind11
     - Topological FCA
     - In-process C, any graph

Contact Process
---------------

The contact process models infection spreading or neural excitation-inhibition.
Two variants are provided:

- **ContactProcessSIR**: Infection/recovery with rate ``mu``
- **ContactProcessEI**: Excitation-inhibition with coupling ``gamma``

ContactProcessEI (Recommended)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Maps to the optimized C1* backends:

.. code-block:: python

   from lrgsglib.graphs import Lattice2D
   from lrgsglib.statsys.ContactProcess import ContactProcessEI

   # Create graph
   lattice = Lattice2D(side1=64, geo='sqr', pflip=0.2, seed=42)
   lattice.flip_random_fract_edges()

   # Initialize contact process
   cp = ContactProcessEI(
       sg=lattice,           # SignedGraph instance
       gamma=1.5,            # Coupling strength (critical point ~1.65)
       steps=5000,           # Number of time steps
       runlang='C1c',        # C backend with snapshot output
       ic='random',          # Initial condition: 'random', 'all_active', 'single'
       rho_init=0.5,         # Initial density for 'random' IC
       seed=42,              # Random seed
   )

   cp.init_contact_dynamics()
   cp.run(verbose=False)

   # Access final state
   active_fraction = cp.s.mean()
   print(f"Final active fraction: {active_fraction:.4f}")

C Backend Options for Contact Process
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table:: Contact Process C Backend Options
   :header-rows: 1
   :widths: 15 85

   * - runlang
     - Description
   * - ``C1a``
     - Final state only
   * - ``C1b``
     - Time series of active fraction
   * - ``C1c``
     - Periodic snapshots of configuration
   * - ``C1d``
     - Cluster statistics
   * - ``C1e``
     - Extended diagnostics

ContactProcessSIR
~~~~~~~~~~~~~~~~~

Alternative parameterization with infection rate:

.. code-block:: python

   from lrgsglib.statsys.ContactProcess import ContactProcessSIR

   cp_sir = ContactProcessSIR(
       sg=lattice,
       mu=0.2,               # Infection rate
       steps=5000,
       runlang='py',         # Python backend for this variant
   )

   cp_sir.init_contact_dynamics()
   cp_sir.run(tqdm_on=True)

Voter Model
-----------

The voter model simulates opinion dynamics where nodes copy neighbors' states.

Basic Usage
~~~~~~~~~~~

.. code-block:: python

   from lrgsglib.graphs import Lattice2D
   from lrgsglib.statsys.VoterModel import VoterModel

   # Create graph
   lattice = Lattice2D(side1=32, geo='sqr', pflip=0.1, seed=42)
   lattice.flip_random_fract_edges()

   # Initialize voter model
   voter = VoterModel(
       sg=lattice,
       steps=10000,
       save_magnetization=True,  # Track magnetization over time
       runlang='C1',             # C backend
   )

   voter.init_voter_dynamics()
   voter.run(verbose=False)

   # Access magnetization trajectory
   print(f"Final magnetization: {voter.magn[-1]:.4f}")
   print(f"Trajectory length: {len(voter.magn)}")

Signed Edge Effects
~~~~~~~~~~~~~~~~~~~

On signed graphs, negative edges cause anti-copying:

.. code-block:: python

   # Node adopts OPPOSITE of neighbor's state across negative edges
   # This creates persistent disagreement patterns

   # Create highly frustrated graph
   frustrated = Lattice2D(side1=32, pflip=0.5, seed=42)
   frustrated.flip_random_fract_edges()

   voter_frust = VoterModel(sg=frustrated, steps=5000, runlang='py')
   voter_frust.init_voter_dynamics()
   voter_frust.run()

Working with Results
--------------------

All dynamics classes store results in attributes after running.

Accessing Observables
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Ising model
   ising.magn       # List of magnetization values
   ising.ene        # List of energy values
   ising.s          # Final spin configuration (numpy array)
   ising.sini       # Initial spin configuration

   # Contact process
   cp.s             # Final state (0/1 for each node)
   cp.s_t           # List of state snapshots (if savedyn=True)

   # Voter model
   voter.magn       # Magnetization trajectory
   voter.s          # Final state

Visualizing Dynamics
~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import matplotlib.pyplot as plt
   import numpy as np

   # Plot Ising magnetization and energy
   fig, axes = plt.subplots(1, 2, figsize=(12, 4))

   axes[0].plot(ising.magn)
   axes[0].set_xlabel('MC Sweep')
   axes[0].set_ylabel('Magnetization')
   axes[0].set_title('Ising Magnetization')

   axes[1].plot(ising.ene)
   axes[1].set_xlabel('MC Sweep')
   axes[1].set_ylabel('Energy')
   axes[1].set_title('Ising Energy')

   plt.tight_layout()
   plt.show()

Energy Calculations
~~~~~~~~~~~~~~~~~~~

Compute energy from spin configurations:

.. code-block:: python

   from lrgsglib.utils.lrg.ising import compute_ising_pairwise_energy

   # Compute energy for current configuration
   energy = ising.calc_full_energy()
   print(f"Current energy: {energy:.4f}")

Batch Simulations
-----------------

For parameter sweeps, use the standalone programs or write loops:

Command-Line Programs
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Ising dynamics on 2D lattice
   python src/L2D_IsingDynamics.py --side 128 --pflip 0.3 --temp 2.0 --steps 10000

   # Contact process scan
   python src/L2D_ContactProcess.py --side 64 --gamma 1.5 --steps 5000

Temperature Sweep Example
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   import numpy as np

   temperatures = np.linspace(1.0, 4.0, 20)
   magnetizations = []

   for T in temperatures:
       ising = IsingDynamics(sg=lattice, T=T, steps=5000, runlang='C1b')
       ising.init_ising_dynamics()
       ising.run(verbose=False)

       # Average magnetization in equilibrium (last 20%)
       eq_magn = np.mean(np.abs(ising.magn[-1000:]))
       magnetizations.append(eq_magn)

   # Plot phase transition
   plt.plot(temperatures, magnetizations, 'o-')
   plt.xlabel('Temperature')
   plt.ylabel('|Magnetization|')
   plt.title('Ising Phase Transition')
   plt.show()

Serialization and Data Management
---------------------------------

Dynamics classes integrate with the library's serialization system:

.. code-block:: python

   # Results are automatically saved to configured paths
   # Check paths:
   print(f"Data directory: {lattice.path_sgdata}")
   print(f"Ising data: {lattice.path_ising}")

   # Output files follow naming conventions:
   # {observable}_{N}_{pflip}_{runid}.bin

Best Practices
--------------

1. **Always use** ``sg=`` **parameter**:

   .. code-block:: python

      # Correct
      IsingDynamics(sg=lattice, T=2.0)

      # Wrong - will fail!
      IsingDynamics(graph=lattice, T=2.0)

2. **Call init method before run**:

   .. code-block:: python

      ising.init_ising_dynamics()  # Required!
      ising.run()

3. **Choose the right backend for your use case**:

   .. code-block:: python

      # GPU (fastest for large systems, any graph)
      runlang='cu_met'

      # Pybind11 (fast, any graph, no file I/O)
      runlang='pb_met'

      # C subprocess (fast, NX graphs only)
      runlang='C1b'

      # Python (slow, for debugging)
      runlang='py'

4. **Set seeds for reproducibility**:

   .. code-block:: python

      lattice = Lattice2D(side1=32, pflip=0.3, seed=42)
      ising = IsingDynamics(sg=lattice, T=2.0, seed=42)

5. **Use appropriate equilibration**:

   .. code-block:: python

      # Discard early samples for equilibration
      equilibrated_magn = ising.magn[len(ising.magn)//5:]

6. **Check backend availability**:

   .. code-block:: python

      # If C backend fails, fall back to Python
      try:
          ising.run(verbose=False)
      except FileNotFoundError:
          print("C backend not found, using Python")
          ising.runlang = 'py'
          ising.run()

API Reference
-------------

.. seealso::

   - :class:`lrgsglib.statsys.IsingDynamics.IsingDynamics` - Ising model
   - :class:`lrgsglib.statsys.ContactProcess.ContactProcessEI` - Contact process
   - :class:`lrgsglib.statsys.ContactProcess.ContactProcessSIR` - SIR variant
   - :class:`lrgsglib.statsys.VoterModel.VoterModel` - Voter model
   - :class:`lrgsglib.statsys.BinDynSys.BinDynSys` - Base class
