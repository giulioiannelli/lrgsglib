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

C Backend Selection
~~~~~~~~~~~~~~~~~~~

The ``runlang`` parameter selects the simulation algorithm:

.. list-table:: Ising C Backend Options
   :header-rows: 1
   :widths: 15 20 65

   * - runlang
     - Algorithm
     - Description
   * - ``C1b``
     - Metropolis
     - Standard Glauber-Metropolis dynamics with E/M output
   * - ``C3b``
     - Simulated Annealing
     - Cooling schedule from high to low temperature
   * - ``C4b``
     - Parallel Tempering
     - Multiple replicas with temperature exchanges

**Naming convention:**
- Number = Algorithm: 1=Metropolis, 3=SA, 4=PT
- Letter = Output: a=final state, b=E/M time series, c=snapshots

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

Python Backend
~~~~~~~~~~~~~~

For debugging or when C backends aren't available:

.. code-block:: python

   ising_py = IsingDynamics(
       sg=lattice,
       T=2.0,
       steps=1000,
       runlang='py',  # Pure Python (slower but always available)
   )

   ising_py.init_ising_dynamics()
   ising_py.run(tqdm_on=True)  # Show progress bar

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

3. **Use C backends for production** (~100x faster):

   .. code-block:: python

      # Fast
      runlang='C1b'

      # Slow (for debugging)
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
