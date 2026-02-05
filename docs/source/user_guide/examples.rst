Complete Examples
=================

This section provides end-to-end examples demonstrating common workflows with
lrgsglib. All examples use fixed seeds for reproducibility.

.. contents:: Contents
   :local:
   :depth: 2

Example 1: Signed Erdos-Renyi Analysis
--------------------------------------

Create a random signed graph and analyze its spectral properties.

.. code-block:: python

   """
   Signed Erdos-Renyi Graph Analysis

   This example creates a random signed network, computes its
   spectrum, and analyzes frustration properties.
   """
   import numpy as np
   import matplotlib.pyplot as plt
   from lrgsglib.graphs.ErdosRenyi import ErdosRenyi
   from lrgsglib.utils.lrg.spectral import get_graph_lspectrum

   # Create signed Erdos-Renyi graph
   # - 200 nodes
   # - Connection probability 0.1
   # - 30% of edges will be negative
   er = ErdosRenyi(nnodes=200, prob=0.1, pflip=0.3, seed=42)
   er.flip_random_fract_edges()

   print(f"Graph created: {er.N} nodes, {er.gr['G'].number_of_edges()} edges")

   # Compute the Laplacian spectrum
   er.compute_laplacian_spectrum()

   # Analyze frustration via negative eigenvalues
   negative_eigvals = er.eigv[er.eigv < -1e-10]
   positive_eigvals = er.eigv[er.eigv > 1e-10]

   print(f"\nSpectral Analysis:")
   print(f"  Total eigenvalues: {len(er.eigv)}")
   print(f"  Negative eigenvalues: {len(negative_eigvals)}")
   print(f"  Zero eigenvalues: {len(er.eigv) - len(negative_eigvals) - len(positive_eigvals)}")
   print(f"  Spectral gap: {er.eigv[1] - er.eigv[0]:.4f}")
   print(f"  Min eigenvalue: {er.eigv.min():.4f}")
   print(f"  Max eigenvalue: {er.eigv.max():.4f}")

   # Plot spectrum
   fig, axes = plt.subplots(1, 2, figsize=(12, 4))

   # Eigenvalue distribution
   axes[0].hist(er.eigv, bins=50, edgecolor='black', alpha=0.7)
   axes[0].axvline(x=0, color='r', linestyle='--', label='Zero')
   axes[0].set_xlabel('Eigenvalue')
   axes[0].set_ylabel('Count')
   axes[0].set_title('Eigenvalue Distribution')
   axes[0].legend()

   # Sorted eigenvalues (spectral staircase)
   axes[1].plot(np.sort(er.eigv), 'b-', linewidth=0.5)
   axes[1].axhline(y=0, color='r', linestyle='--', alpha=0.5)
   axes[1].set_xlabel('Index')
   axes[1].set_ylabel('Eigenvalue')
   axes[1].set_title('Sorted Eigenvalues')

   plt.tight_layout()
   plt.savefig('erdos_renyi_spectrum.png', dpi=150)
   plt.show()

Example 2: Contact Process on 2D Lattice
----------------------------------------

Simulate excitation-inhibition dynamics near the critical point.

.. code-block:: python

   """
   Contact Process Simulation

   This example runs a contact process on a 2D lattice
   with signed edges and analyzes the critical behavior.
   """
   import numpy as np
   import matplotlib.pyplot as plt
   from lrgsglib.graphs.Lattice2D import Lattice2D
   from lrgsglib.statsys.ContactProcess import ContactProcessEI

   # Create a 64x64 square lattice with 20% frustrated edges
   lattice = Lattice2D(side1=64, geo='sqr', pflip=0.2, seed=42)
   lattice.flip_random_fract_edges()

   print(f"Lattice: {lattice.N} nodes, pflip={lattice.pflip}")

   # Scan gamma values around critical point
   gamma_values = np.linspace(1.0, 2.5, 15)
   final_densities = []

   for gamma in gamma_values:
       cp = ContactProcessEI(
           sg=lattice,
           gamma=gamma,
           steps=2000,
           runlang='py',  # Use Python backend (C1c for production)
           ic='random',
           rho_init=0.5,
           seed=42,
       )
       cp.init_contact_dynamics()
       cp.run(verbose=False)

       # Final active fraction
       density = cp.s.mean()
       final_densities.append(density)
       print(f"gamma={gamma:.2f}, final density={density:.3f}")

   # Plot phase diagram
   plt.figure(figsize=(8, 5))
   plt.plot(gamma_values, final_densities, 'o-', markersize=8)
   plt.axvline(x=1.65, color='r', linestyle='--', label='Critical point (~1.65)')
   plt.xlabel('Coupling strength (gamma)')
   plt.ylabel('Final active fraction')
   plt.title('Contact Process Phase Diagram')
   plt.legend()
   plt.grid(True, alpha=0.3)
   plt.savefig('contact_process_phase.png', dpi=150)
   plt.show()

Example 3: Laplacian Spectrum and Entropy
-----------------------------------------

Compute entropy observables from the full Laplacian spectrum.

.. code-block:: python

   """
   Spectral Entropy Analysis

   This example computes the von Neumann entropy profile
   and estimates the spectral dimension.
   """
   import numpy as np
   import matplotlib.pyplot as plt
   from lrgsglib.graphs.MultiplicativeCascade import MultiplicativeCascadeGraph
   from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
   from lrgsglib.utils.lrg.infocomm import compute_entropy_observables_from_eigenvalues

   # Create a multiplicative cascade graph (hierarchical structure)
   mc = MultiplicativeCascadeGraph(
       p1=0.8, p2=0.6, p3=0.6, p4=0.8,
       fraction=0.4,
       iterations=6,
       seed=42,
   )

   print(f"Graph: {mc.N} nodes")

   # Compute full spectrum (required for entropy)
   L, eigenvalues = get_graph_lspectrum(mc.gr['G'], library='numpy')
   eigenvalues = np.sort(np.real(eigenvalues))  # Ensure sorted, real

   print(f"Spectrum computed: {len(eigenvalues)} eigenvalues")
   print(f"Spectral range: [{eigenvalues.min():.4f}, {eigenvalues.max():.4f}]")

   # Compute entropy observables
   entropy, specific_heat, variance, time_grid = compute_entropy_observables_from_eigenvalues(
       eigenvalues=eigenvalues,
       num_nodes=mc.N,
       steps=300,
       t1=-2,
       t2=5,
       entropy_norm='complement',  # Returns 1 - S/log(N)
       specific_heat_scale='logN',
   )

   # Estimate spectral dimension from tail of specific heat
   tail_fraction = 0.1
   tail_len = int(tail_fraction * len(specific_heat))
   ds_estimate = 2.0 * np.mean(specific_heat[-tail_len:])

   print(f"\nResults:")
   print(f"  Spectral dimension estimate: {ds_estimate:.3f}")

   # Plot results
   fig, axes = plt.subplots(1, 3, figsize=(15, 4))

   # Entropy profile
   axes[0].semilogx(time_grid, entropy)
   axes[0].set_xlabel(r'$\tau$')
   axes[0].set_ylabel(r'$1 - S(\tau)/\log N$')
   axes[0].set_title('Entropy Profile')
   axes[0].grid(True, alpha=0.3)

   # Specific heat (spectral dimension proxy)
   axes[1].semilogx(time_grid, specific_heat)
   axes[1].axhline(y=ds_estimate/2, color='r', linestyle='--',
                   label=f'$d_s/2 = {ds_estimate/2:.2f}$')
   axes[1].set_xlabel(r'$\tau$')
   axes[1].set_ylabel(r'$C(\tau)$')
   axes[1].set_title('Specific Heat')
   axes[1].legend()
   axes[1].grid(True, alpha=0.3)

   # Variance
   axes[2].loglog(time_grid, variance)
   axes[2].set_xlabel(r'$\tau$')
   axes[2].set_ylabel('Variance')
   axes[2].set_title('Eigenvalue Variance')
   axes[2].grid(True, alpha=0.3)

   plt.tight_layout()
   plt.savefig('entropy_analysis.png', dpi=150)
   plt.show()

Example 4: Ising Dynamics and Phase Transition
----------------------------------------------

Study the ferromagnetic-paramagnetic phase transition on a frustrated lattice.

.. code-block:: python

   """
   Ising Model Phase Transition

   This example explores the phase transition on a 2D lattice
   with variable frustration.
   """
   import numpy as np
   import matplotlib.pyplot as plt
   from lrgsglib.graphs.Lattice2D import Lattice2D
   from lrgsglib.statsys.IsingDynamics import IsingDynamics

   # Create 32x32 lattice with moderate frustration
   lattice = Lattice2D(side1=32, geo='sqr', pflip=0.15, seed=42)
   lattice.flip_random_fract_edges()

   print(f"Lattice: {lattice.N} nodes, pflip={lattice.pflip}")

   # Temperature sweep
   temperatures = np.linspace(1.0, 4.0, 20)
   magnetizations = []
   energies = []

   for T in temperatures:
       ising = IsingDynamics(
           sg=lattice,
           T=T,
           steps=5000,
           thrmSTEP=1000,  # Thermalization steps
           runlang='py',   # Use 'C1b' for production
       )
       ising.init_ising_dynamics()
       ising.run(verbose=False)

       # Equilibrium averages (last 20% of samples)
       eq_start = len(ising.magn) * 4 // 5
       avg_magn = np.mean(np.abs(ising.magn[eq_start:]))
       avg_ene = np.mean(ising.ene[eq_start:]) if ising.ene else 0

       magnetizations.append(avg_magn)
       energies.append(avg_ene)
       print(f"T={T:.2f}: |M|={avg_magn:.3f}, E={avg_ene:.1f}")

   # Plot phase transition
   fig, axes = plt.subplots(1, 2, figsize=(12, 4))

   # Magnetization vs Temperature
   axes[0].plot(temperatures, magnetizations, 'o-', markersize=6)
   axes[0].set_xlabel('Temperature')
   axes[0].set_ylabel('|Magnetization|')
   axes[0].set_title('Order Parameter')
   axes[0].grid(True, alpha=0.3)

   # Energy vs Temperature
   axes[1].plot(temperatures, energies, 's-', color='red', markersize=6)
   axes[1].set_xlabel('Temperature')
   axes[1].set_ylabel('Energy')
   axes[1].set_title('Internal Energy')
   axes[1].grid(True, alpha=0.3)

   plt.tight_layout()
   plt.savefig('ising_phase_transition.png', dpi=150)
   plt.show()

   # Time series at critical temperature
   T_c = 2.27  # Approximate critical temperature
   ising_tc = IsingDynamics(sg=lattice, T=T_c, steps=10000, runlang='py')
   ising_tc.init_ising_dynamics()
   ising_tc.run(verbose=False)

   plt.figure(figsize=(10, 3))
   plt.plot(ising_tc.magn, linewidth=0.5)
   plt.xlabel('MC Sweep')
   plt.ylabel('Magnetization')
   plt.title(f'Magnetization Time Series at T={T_c}')
   plt.savefig('ising_timeseries.png', dpi=150)
   plt.show()

Example 5: Comparing Python and C Backends
------------------------------------------

Demonstrate the performance difference between Python and C backends.

.. code-block:: python

   """
   Backend Performance Comparison

   This example compares the speed of Python vs C backends
   for Ising dynamics simulation.
   """
   import time
   import numpy as np
   from lrgsglib.graphs.Lattice2D import Lattice2D
   from lrgsglib.statsys.IsingDynamics import IsingDynamics

   # Create lattice
   lattice = Lattice2D(side1=32, geo='sqr', pflip=0.2, seed=42)
   lattice.flip_random_fract_edges()

   print(f"Benchmarking on {lattice.N}-node lattice")
   print("=" * 50)

   # Python backend
   print("\nPython backend (runlang='py'):")
   ising_py = IsingDynamics(sg=lattice, T=2.0, steps=1000, runlang='py')
   ising_py.init_ising_dynamics()

   start = time.time()
   ising_py.run(verbose=False)
   py_time = time.time() - start
   print(f"  Time: {py_time:.3f}s")
   print(f"  Final magnetization: {ising_py.magn[-1]:.4f}")

   # C backend (if available)
   print("\nC backend (runlang='C1b'):")
   try:
       ising_c = IsingDynamics(sg=lattice, T=2.0, steps=1000, runlang='C1b')
       ising_c.init_ising_dynamics()

       start = time.time()
       ising_c.run(verbose=False)
       c_time = time.time() - start
       print(f"  Time: {c_time:.3f}s")
       print(f"  Final magnetization: {ising_c.magn[-1]:.4f}")
       print(f"\nSpeedup: {py_time/c_time:.1f}x")
   except FileNotFoundError:
       print("  C backend not available (not compiled)")
       print("  Build C extensions with: cd lrgsglib && make all")

   # Spectral computation comparison
   print("\n" + "=" * 50)
   print("Spectral Computation:")

   from lrgsglib.utils.lrg.spectral import get_graph_lspectrum

   for backend in ['numpy', 'scipy']:
       start = time.time()
       L, w = get_graph_lspectrum(lattice.gr['G'], library=backend)
       elapsed = time.time() - start
       print(f"  {backend}: {elapsed:.3f}s")

Example 6: 3D Lattice Analysis
------------------------------

Analyze spectral properties of a 3D cubic lattice.

.. code-block:: python

   """
   3D Lattice Spectral Analysis

   This example creates a 3D cubic lattice and analyzes
   its spectral properties under frustration.
   """
   import numpy as np
   import matplotlib.pyplot as plt
   from lrgsglib.graphs.Lattice3D import Lattice3D
   from lrgsglib.utils.lrg.infocomm import compute_entropy_observables_from_eigenvalues

   # Create 10x10x10 cubic lattice with signed edges
   lattice3d = Lattice3D(
       side1=10,
       geo='cub',     # cubic geometry
       pflip=0.1,     # 10% frustrated edges
       seed=42,
   )
   lattice3d.flip_random_fract_edges()

   print(f"3D Lattice: {lattice3d.N} nodes ({lattice3d.side1}^3)")
   print(f"Edges: {lattice3d.gr['G'].number_of_edges()}")

   # Compute spectrum
   lattice3d.compute_laplacian_spectrum()

   # Compare signed vs unsigned spectrum
   n_negative = np.sum(lattice3d.eigv < -1e-10)
   print(f"\nSpectral Analysis:")
   print(f"  Negative eigenvalues: {n_negative}")
   print(f"  Spectral range: [{lattice3d.eigv.min():.4f}, {lattice3d.eigv.max():.4f}]")

   # Compute entropy
   entropy, heat, var, tau = compute_entropy_observables_from_eigenvalues(
       eigenvalues=lattice3d.eigv,
       num_nodes=lattice3d.N,
       steps=200,
   )

   # Spectral dimension
   tail_len = int(0.1 * len(heat))
   ds = 2.0 * np.mean(heat[-tail_len:])
   print(f"  Spectral dimension: {ds:.2f} (expected ~3 for 3D lattice)")

   # Plot eigenvalue density for different pflip values
   fig, axes = plt.subplots(1, 2, figsize=(12, 4))

   # Eigenvalue histogram
   axes[0].hist(lattice3d.eigv, bins=50, density=True, alpha=0.7, edgecolor='black')
   axes[0].axvline(x=0, color='r', linestyle='--')
   axes[0].set_xlabel('Eigenvalue')
   axes[0].set_ylabel('Density')
   axes[0].set_title(f'3D Lattice Spectrum (pflip={lattice3d.pflip})')

   # Specific heat
   axes[1].semilogx(tau, heat)
   axes[1].axhline(y=ds/2, color='r', linestyle='--', label=f'$d_s/2={ds/2:.2f}$')
   axes[1].set_xlabel(r'$\tau$')
   axes[1].set_ylabel(r'$C(\tau)$')
   axes[1].set_title('Specific Heat')
   axes[1].legend()
   axes[1].grid(True, alpha=0.3)

   plt.tight_layout()
   plt.savefig('lattice3d_analysis.png', dpi=150)
   plt.show()

   # Compare different frustration levels
   print("\nFrustration comparison:")
   for pflip in [0.0, 0.1, 0.3, 0.5]:
       lat = Lattice3D(side1=8, geo='cub', pflip=pflip, seed=42)
       lat.flip_random_fract_edges()
       lat.compute_laplacian_spectrum()
       n_neg = np.sum(lat.eigv < -1e-10)
       print(f"  pflip={pflip}: {n_neg} negative eigenvalues ({100*n_neg/lat.N:.1f}%)")

Example 7: Custom Analysis Pipeline
-----------------------------------

Build a complete analysis pipeline for signed graph research.

.. code-block:: python

   """
   Custom Analysis Pipeline

   This example demonstrates a complete workflow for
   analyzing signed networks from creation to visualization.
   """
   import numpy as np
   import matplotlib.pyplot as plt
   from lrgsglib.graphs.Lattice2D import Lattice2D
   from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
   from lrgsglib.utils.lrg.infocomm import (
       compute_entropy_observables_from_eigenvalues,
       compute_renyi_observables_from_eigenvalues,
   )
   from lrgsglib.statsys.IsingDynamics import IsingDynamics

   # =========================================================
   # Step 1: Create and characterize the graph
   # =========================================================
   print("Step 1: Graph Creation")
   print("-" * 40)

   lattice = Lattice2D(side1=32, geo='sqr', pflip=0.25, seed=42)
   lattice.flip_random_fract_edges()

   # Count positive and negative edges
   G = lattice.gr['G']
   n_positive = sum(1 for u, v in G.edges() if G[u][v].get('weight', 1) > 0)
   n_negative = G.number_of_edges() - n_positive

   print(f"Nodes: {lattice.N}")
   print(f"Edges: {G.number_of_edges()} (positive: {n_positive}, negative: {n_negative})")
   print(f"Frustration parameter: pflip={lattice.pflip}")

   # =========================================================
   # Step 2: Spectral analysis
   # =========================================================
   print("\nStep 2: Spectral Analysis")
   print("-" * 40)

   L, eigenvalues = get_graph_lspectrum(G, library='numpy')
   eigenvalues = np.sort(np.real(eigenvalues))

   spectral_gap = eigenvalues[1] - eigenvalues[0]
   n_negative_eig = np.sum(eigenvalues < -1e-10)

   print(f"Spectral gap: {spectral_gap:.6f}")
   print(f"Negative eigenvalues: {n_negative_eig}")
   print(f"Min eigenvalue: {eigenvalues.min():.4f}")

   # =========================================================
   # Step 3: Entropy and information measures
   # =========================================================
   print("\nStep 3: Entropy Analysis")
   print("-" * 40)

   # Von Neumann entropy
   entropy, heat, var, tau = compute_entropy_observables_from_eigenvalues(
       eigenvalues=eigenvalues,
       num_nodes=lattice.N,
       steps=200,
   )

   # Spectral dimension
   tail_len = int(0.1 * len(heat))
   ds = 2.0 * np.mean(heat[-tail_len:])

   # Renyi entropy (q=2)
   renyi_results = compute_renyi_observables_from_eigenvalues(
       eigenvalues=eigenvalues,
       num_nodes=lattice.N,
       q=2.0,
       steps=200,
   )

   print(f"Spectral dimension: {ds:.3f}")
   print(f"Renyi-2 ds estimate: {renyi_results['ds_estimate']:.3f}")

   # =========================================================
   # Step 4: Dynamical simulations
   # =========================================================
   print("\nStep 4: Ising Dynamics")
   print("-" * 40)

   # Run at critical temperature
   T_c = 2.27
   ising = IsingDynamics(sg=lattice, T=T_c, steps=5000, runlang='py')
   ising.init_ising_dynamics()
   ising.run(verbose=False)

   # Compute observables
   eq_magn = ising.magn[len(ising.magn)//2:]
   avg_magn = np.mean(np.abs(eq_magn))
   std_magn = np.std(eq_magn)

   print(f"Temperature: {T_c}")
   print(f"Average |M|: {avg_magn:.4f}")
   print(f"M fluctuation: {std_magn:.4f}")

   # =========================================================
   # Step 5: Comprehensive visualization
   # =========================================================
   print("\nStep 5: Visualization")
   print("-" * 40)

   fig = plt.figure(figsize=(15, 10))

   # Eigenvalue spectrum
   ax1 = fig.add_subplot(2, 3, 1)
   ax1.hist(eigenvalues, bins=50, edgecolor='black', alpha=0.7)
   ax1.axvline(x=0, color='r', linestyle='--')
   ax1.set_xlabel('Eigenvalue')
   ax1.set_ylabel('Count')
   ax1.set_title('Laplacian Spectrum')

   # Entropy profile
   ax2 = fig.add_subplot(2, 3, 2)
   ax2.semilogx(tau, entropy)
   ax2.set_xlabel(r'$\tau$')
   ax2.set_ylabel(r'$1 - S(\tau)/\log N$')
   ax2.set_title('Entropy Profile')
   ax2.grid(True, alpha=0.3)

   # Specific heat
   ax3 = fig.add_subplot(2, 3, 3)
   ax3.semilogx(tau, heat)
   ax3.axhline(y=ds/2, color='r', linestyle='--', label=f'$d_s/2={ds/2:.2f}$')
   ax3.set_xlabel(r'$\tau$')
   ax3.set_ylabel(r'$C(\tau)$')
   ax3.set_title('Specific Heat')
   ax3.legend()
   ax3.grid(True, alpha=0.3)

   # Ising magnetization
   ax4 = fig.add_subplot(2, 3, 4)
   ax4.plot(ising.magn, linewidth=0.5)
   ax4.set_xlabel('MC Sweep')
   ax4.set_ylabel('Magnetization')
   ax4.set_title(f'Ising at T={T_c}')

   # Magnetization histogram
   ax5 = fig.add_subplot(2, 3, 5)
   ax5.hist(eq_magn, bins=50, edgecolor='black', alpha=0.7)
   ax5.set_xlabel('Magnetization')
   ax5.set_ylabel('Count')
   ax5.set_title('Equilibrium Distribution')

   # Renyi vs Von Neumann entropy
   ax6 = fig.add_subplot(2, 3, 6)
   ax6.semilogx(tau, entropy, label='Von Neumann')
   ax6.semilogx(renyi_results['tau'], renyi_results['renyi_entropy'],
                label='Renyi (q=2)', linestyle='--')
   ax6.set_xlabel(r'$\tau$')
   ax6.set_ylabel('Entropy')
   ax6.set_title('Entropy Comparison')
   ax6.legend()
   ax6.grid(True, alpha=0.3)

   plt.tight_layout()
   plt.savefig('complete_analysis.png', dpi=150)
   plt.show()

   # =========================================================
   # Step 6: Summary report
   # =========================================================
   print("\nSummary Report")
   print("=" * 50)
   print(f"Graph: {lattice.N} nodes, {G.number_of_edges()} edges")
   print(f"Frustration: pflip={lattice.pflip}, {n_negative_eig} negative eigenvalues")
   print(f"Spectral dimension: {ds:.3f}")
   print(f"Critical Ising magnetization: {avg_magn:.4f} +/- {std_magn:.4f}")
   print("=" * 50)

Running the Examples
--------------------

All examples can be run as standalone Python scripts. To run them:

.. code-block:: bash

   # Ensure lrgsglib is installed
   pip install -e lrgsglib

   # Run an example
   python example_script.py

For best performance, ensure C backends are compiled:

.. code-block:: bash

   cd lrgsglib && make all

Common Issues
-------------

1. **ImportError**: Ensure lrgsglib is installed in editable mode
2. **FileNotFoundError for C backend**: Compile C extensions with ``make all``
3. **Slow execution**: Use C backends (``runlang='C1b'``) instead of Python
4. **Memory errors**: Reduce graph size or use sparse methods

.. seealso::

   - :doc:`spectral` - Detailed spectral analysis guide
   - :doc:`dynamics` - Statistical physics simulations
   - :doc:`graphs` - Graph creation and manipulation
   - :doc:`plotting` - Visualization techniques
