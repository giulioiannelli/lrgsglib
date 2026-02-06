Statistical Physics on Networks
================================

This section describes the statistical physics models implemented in lrgsglib:
Ising model, contact process, and voter model.

Ising Model
-----------

The **Ising model** describes interacting spins on a network.

Hamiltonian
~~~~~~~~~~~

For a signed graph with spins :math:`s_i \in \{-1, +1\}`:

.. math::

   H = -\sum_{(i,j) \in E} J_{ij} s_i s_j

where :math:`J_{ij}` is the coupling strength:

- :math:`J_{ij} > 0`: Ferromagnetic (prefer alignment)
- :math:`J_{ij} < 0`: Antiferromagnetic (prefer anti-alignment)

For signed graphs with unit couplings: :math:`J_{ij} = A_{ij}` (signed adjacency).

Boltzmann Distribution
~~~~~~~~~~~~~~~~~~~~~~

At temperature T, the probability of configuration s:

.. math::

   P(s) = \frac{e^{-\beta H(s)}}{Z}

where:

- :math:`\beta = 1/T` is inverse temperature
- :math:`Z = \sum_s e^{-\beta H(s)}` is the partition function

Observables
~~~~~~~~~~~

**Magnetization:**

.. math::

   m = \frac{1}{N} \sum_i s_i

**Energy per spin:**

.. math::

   e = \frac{1}{N} \langle H \rangle

**Susceptibility:**

.. math::

   \chi = N(\langle m^2 \rangle - \langle m \rangle^2)

Phase Transition
~~~~~~~~~~~~~~~~

In the thermodynamic limit:

- **High T** (β → 0): Disordered, m ≈ 0
- **Low T** (β → ∞): Ordered, |m| > 0
- **Critical T = Tc**: Phase transition

For 2D square lattice without frustration: :math:`T_c = 2/\ln(1+\sqrt{2}) \approx 2.269`

Frustration lowers Tc and can destroy long-range order.

Monte Carlo Methods
-------------------

Metropolis Algorithm
~~~~~~~~~~~~~~~~~~~~

Single-spin update at temperature T:

1. Select random spin i
2. Compute energy change: :math:`\Delta E = 2 s_i \sum_j J_{ij} s_j`
3. Accept flip with probability :math:`\min(1, e^{-\beta \Delta E})`

Equilibration requires O(N) single-spin updates = 1 Monte Carlo sweep.

.. code-block:: python

   from lrgsglib.statsys.IsingDynamics import IsingDynamics
   from lrgsglib.graphs.Lattice2D import Lattice2D

   # Create frustrated lattice
   lattice = Lattice2D(side=32, pflip=0.3, seed=42)
   lattice.flip_random_fract_edges()

   # Metropolis simulation
   ising = IsingDynamics(sg=lattice, T=2.0, steps=10000, runlang='C1b')
   ising.init_ising_dynamics()
   ising.run()

   print(f"Final magnetization: {ising.magn[-1]:.4f}")

Simulated Annealing
~~~~~~~~~~~~~~~~~~~

For finding ground states, gradually lower temperature:

.. math::

   T(t) = T_0 \cdot \alpha^t

where α < 1 is the cooling rate.

.. code-block:: python

   # Simulated annealing
   ising = IsingDynamics(
       sg=lattice,
       T=5.0,           # Starting temperature
       Tfin=0.01,       # Final temperature
       steps=50000,
       runlang='C3b'    # SA backend
   )
   ising.init_ising_dynamics()
   ising.run()

Parallel Tempering
~~~~~~~~~~~~~~~~~~

Run multiple replicas at different temperatures and exchange configurations:

1. Run each replica with Metropolis
2. Attempt swap between adjacent temperatures
3. Accept swap with probability based on detailed balance

.. code-block:: python

   # Parallel tempering
   ising = IsingDynamics(
       sg=lattice,
       T=1.0,           # Lowest temperature
       Tmax=4.0,        # Highest temperature
       n_replicas=8,
       steps=10000,
       runlang='C4b'    # PT backend
   )
   ising.init_ising_dynamics()
   ising.run()

Contact Process
---------------

The **contact process** models epidemic spreading on networks.

Dynamics
~~~~~~~~

Each node is either:

- **Infected** (I): State 1
- **Susceptible** (S): State 0

Transitions:

- **Infection**: S → I at rate λ·(number of infected neighbors)
- **Recovery**: I → S at rate μ (often set to 1)

The control parameter is :math:`\gamma = \lambda/\mu` (infection rate).

Critical Behavior
~~~~~~~~~~~~~~~~~

The contact process exhibits a **phase transition**:

- :math:`\gamma < \gamma_c`: Absorbing phase (epidemic dies out)
- :math:`\gamma > \gamma_c`: Active phase (endemic state)

Near criticality:

.. math::

   \rho_{\infty} \sim (\gamma - \gamma_c)^{\beta}

where ρ∞ is the steady-state density of infected nodes.

Implementation
~~~~~~~~~~~~~~

lrgsglib provides two variants:

**EI Model** (Exposed → Infected):

.. code-block:: python

   from lrgsglib.statsys.ContactProcess import ContactProcessEI

   cp = ContactProcessEI(
       sg=lattice,
       gamma=1.5,      # Infection rate
       steps=5000,
       runlang='C1c'
   )
   cp.init_contact_dynamics()
   cp.run()

   # Density time series
   rho_t = cp.rho_t

**SIR Model** (Susceptible → Infected → Recovered):

.. code-block:: python

   from lrgsglib.statsys.ContactProcess import ContactProcessSIR

   cp = ContactProcessSIR(
       sg=lattice,
       gamma=1.5,      # Infection rate
       steps=5000,
       runlang='C1c'
   )
   cp.init_contact_dynamics()
   cp.run()

Voter Model
-----------

The **voter model** describes opinion dynamics where agents adopt neighbors' opinions.

Dynamics
~~~~~~~~

Each node has opinion :math:`s_i \in \{-1, +1\}`.

Update rule (asynchronous):

1. Select random node i
2. Select random neighbor j
3. If edge (i,j) is positive: :math:`s_i \leftarrow s_j`
4. If edge (i,j) is negative: :math:`s_i \leftarrow -s_j`

On unsigned graphs, the system reaches consensus. On signed graphs with
frustration, dynamics may not reach consensus.

Observables
~~~~~~~~~~~

**Opinion density:**

.. math::

   \rho = \frac{1}{N} \sum_i \frac{1 + s_i}{2}

**Magnetization:**

.. math::

   m = \frac{1}{N} \sum_i s_i

Implementation
~~~~~~~~~~~~~~

.. code-block:: python

   from lrgsglib.statsys.VoterModel import VoterModel

   voter = VoterModel(
       sg=lattice,
       steps=10000,
       runlang='C1'
   )
   voter.init_voter_dynamics()
   voter.run()

   # Magnetization time series
   magn_t = voter.magn

Connection to Signed Graphs
---------------------------

Frustrated Ground States
~~~~~~~~~~~~~~~~~~~~~~~~

For the Ising model on signed graphs:

- **Balanced graph**: Unique ground state (up to global flip)
- **Frustrated graph**: Multiple degenerate ground states

The ground state energy:

.. math::

   E_0 = -\sum_{(i,j)} |J_{ij}| + 2 \phi(G)

where φ(G) is the frustration index (number of frustrated edges).

Epidemic Threshold
~~~~~~~~~~~~~~~~~~

The critical infection rate γc depends on network structure:

- Regular lattices: Universal value depending on dimension
- Random graphs: :math:`\gamma_c \sim 1/\langle k \rangle` (mean degree)
- Signed graphs: Frustration can modify the threshold

Spectral Connection
~~~~~~~~~~~~~~~~~~~

For diffusion-like processes, the relaxation time:

.. math::

   \tau \sim \frac{1}{\lambda_2}

where λ₂ is the spectral gap. Frustrated graphs with small λ₂
have slow dynamics.

Finite-Size Effects
-------------------

In finite systems:

- True phase transitions become crossovers
- Susceptibility peaks near Tc
- Correlation length is bounded by system size

Scaling analysis:

.. math::

   m(L, T) = L^{-\beta/\nu} f[(T - T_c) L^{1/\nu}]

where L is system size and ν, β are critical exponents.

C Backend Performance
---------------------

lrgsglib provides C implementations (~100x faster than Python):

.. list-table::
   :header-rows: 1

   * - Model
     - Backend
     - Algorithm
   * - Ising
     - ``C1b``
     - Optimized Metropolis
   * - Ising
     - ``C3b``
     - Simulated Annealing
   * - Ising
     - ``C4b``
     - Parallel Tempering
   * - Contact Process
     - ``C1c``
     - Gillespie algorithm
   * - Voter
     - ``C1``
     - Asynchronous update

Example: Phase Diagram
----------------------

Mapping the Ising phase transition:

.. code-block:: python

   import numpy as np
   from lrgsglib.statsys.IsingDynamics import IsingDynamics
   from lrgsglib.graphs.Lattice2D import Lattice2D

   temperatures = np.linspace(1.0, 4.0, 20)
   magnetizations = []

   for T in temperatures:
       lattice = Lattice2D(side=32, pflip=0.0, seed=42)

       ising = IsingDynamics(sg=lattice, T=T, steps=5000, runlang='C1b')
       ising.init_ising_dynamics()
       ising.run()

       # Average over last half of trajectory
       m_avg = np.mean(np.abs(ising.magn[2500:]))
       magnetizations.append(m_avg)

   # Plot m vs T to see phase transition

See Also
--------

- :doc:`lrg_theory` - Spectral analysis framework
- :doc:`signed_graphs` - Signed graph fundamentals
- :doc:`/user_guide/dynamics` - Practical simulation guide
