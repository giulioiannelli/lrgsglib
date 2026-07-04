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
   from lrgsglib.graphs import ErdosRenyi
   from lrgsglib.utils.lrg.spectral import get_graph_lspectrum

   # Create signed Erdos-Renyi graph
   # - 200 nodes
   # - Connection probability 0.1
   # - 30% of edges will be negative
   er = ErdosRenyi(n=200, p=0.1, pflip=0.3, seed=42)
   er.flip_random_fract_edges()

   print(f"Graph created: {er.N} nodes, {er.gr['G'].number_of_edges()} edges")

   # Compute the Laplacian spectrum
   er.compute_laplacian_spectrum()

   # Frustration signal: the smallest eigenvalue of the (positive-semidefinite)
   # signed Laplacian. lambda_min == 0 -> balanced; lambda_min > 0 -> frustrated.
   eigvals = np.sort(er.eigv)

   print(f"\nSpectral Analysis:")
   print(f"  Total eigenvalues: {len(eigvals)}")
   print(f"  Smallest eigenvalue (frustration signal): {eigvals[0]:.4f}")
   print(f"  Balanced: {np.isclose(eigvals[0], 0.0, atol=1e-8)}")
   print(f"  Spectral gap: {eigvals[1] - eigvals[0]:.4f}")
   print(f"  Min eigenvalue: {eigvals.min():.4f}")
   print(f"  Max eigenvalue: {eigvals.max():.4f}")

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
   from lrgsglib.graphs import Lattice2D
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
   from lrgsglib.graphs.nx import MultiplicativeCascadeGraphNX as MultiplicativeCascadeGraph
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
   from lrgsglib.graphs import Lattice2D
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
   from lrgsglib.graphs import Lattice2D
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
   from lrgsglib.graphs import Lattice3D
   from lrgsglib.utils.lrg.infocomm import compute_entropy_observables_from_eigenvalues

   # Create 10x10x10 cubic lattice with signed edges
   lattice3d = Lattice3D(
       dim=10,
       geo='sc',      # simple-cubic geometry
       pflip=0.1,     # 10% frustrated edges
       seed=42,
   )
   lattice3d.flip_random_fract_edges()

   print(f"3D Lattice: {lattice3d.N} nodes ({lattice3d.dim[0]}^3)")
   print(f"Edges: {lattice3d.gr['G'].number_of_edges()}")

   # Compute spectrum
   lattice3d.compute_laplacian_spectrum()

   # Frustration signal: smallest eigenvalue of the (PSD) signed Laplacian
   # (0 if balanced, > 0 if frustrated).
   eigvals = np.sort(lattice3d.eigv)
   print(f"\nSpectral Analysis:")
   print(f"  Smallest eigenvalue (frustration signal): {eigvals[0]:.4f}")
   print(f"  Spectral range: [{eigvals.min():.4f}, {eigvals.max():.4f}]")

   # Compute entropy
   entropy, heat, var, tau = compute_entropy_observables_from_eigenvalues(
       eigenvalues=lattice3d.eigv,
       num_nodes=lattice3d.N,
       steps=200,
   )

   # Specific-heat peak: approaches d_s/2 in the large-N scaling regime, so
   # 2 * C_max is a finite-size proxy for the spectral dimension (~3 for a 3D
   # lattice only in the thermodynamic limit; a small lattice overshoots).
   c_peak = float(np.max(heat))
   print(f"  Specific-heat peak C_max = {c_peak:.2f}  (~ d_s/2 in the scaling limit)")

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
   axes[1].axhline(y=c_peak, color='r', linestyle='--', label=f'$C_{{max}}={c_peak:.2f}$')
   axes[1].set_xlabel(r'$\tau$')
   axes[1].set_ylabel(r'$C(\tau)$')
   axes[1].set_title('Specific Heat')
   axes[1].legend()
   axes[1].grid(True, alpha=0.3)

   plt.tight_layout()
   plt.savefig('lattice3d_analysis.png', dpi=150)
   plt.show()

   # Compare different frustration levels via the smallest eigenvalue
   # (0 for the balanced pflip=0 case, growing as frustration increases).
   print("\nFrustration comparison (smallest signed-Laplacian eigenvalue):")
   for pflip in [0.0, 0.1, 0.3, 0.5]:
       lat = Lattice3D(dim=8, geo='sc', pflip=pflip, seed=42)
       lat.flip_random_fract_edges()
       lat.compute_laplacian_spectrum()
       lam_min = float(np.sort(lat.eigv)[0])
       print(f"  pflip={pflip}: lambda_min = {lam_min:.4f}")

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
   from lrgsglib.graphs import Lattice2D
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
   lambda_min = eigenvalues[0]   # frustration signal (0 balanced, >0 frustrated)

   print(f"Spectral gap: {spectral_gap:.6f}")
   print(f"Smallest eigenvalue (frustration signal): {lambda_min:.4f}")
   print(f"Max eigenvalue: {eigenvalues.max():.4f}")

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
   print(f"Frustration: pflip={lattice.pflip}, lambda_min = {lambda_min:.4f}")
   print(f"Spectral dimension: {ds:.3f}")
   print(f"Critical Ising magnetization: {avg_magn:.4f} +/- {std_magn:.4f}")
   print("=" * 50)

Example 8: GraphOfGraphs - Hierarchical Networks
-------------------------------------------------

Build hierarchical "graph of graphs" structures with arbitrary base and fiber components.

.. code-block:: python

   """
   GraphOfGraphs: Hierarchical Network Construction

   This example demonstrates creating hierarchical structures
   where a "fiber" graph is attached to each node of a "base" graph.
   """
   import numpy as np
   import matplotlib.pyplot as plt
   from lrgsglib.graphs import GraphOfGraphs, Lattice2D, ErdosRenyi

   # =========================================================
   # Example 1: Basic GraphOfGraphs creation
   # =========================================================
   print("Example 1: Basic Creation")
   print("-" * 40)

   # Create a 2D lattice base with small lattice fibers
   gog = GraphOfGraphs(
       base_graph_type='Lattice2D',
       base_params={'side1': 5, 'geo': 'sqr'},      # 5x5 = 25 base nodes
       fiber_graph_type='Lattice2D',
       fiber_params={'side1': 4, 'geo': 'sqr'},     # 4x4 = 16 fiber nodes each
       anchor_policy='first',                        # Connect fiber node 0 to base
       pflip=0.1,                                    # 10% frustrated edges
       seed=42,
       engine='nx'                                   # Use NetworkX backend
   )

   print(f"Base nodes: {gog.N_base}")
   print(f"Fiber nodes per base: {gog.N_fiber}")
   print(f"Total nodes: {gog.N}")
   print(f"Expected total: {gog.N_base} + {gog.N_base} * {gog.N_fiber} = {gog.N_base + gog.N_base * gog.N_fiber}")

   # =========================================================
   # Example 2: Different anchor policies
   # =========================================================
   print("\nExample 2: Anchor Policies")
   print("-" * 40)

   policies = ['first', 'center', 'last', 'random']
   for policy in policies:
       gog_p = GraphOfGraphs(
           base_graph_type='Lattice2D',
           base_params={'side1': 3},
           fiber_graph_type='Lattice2D',
           fiber_params={'side1': 3},
           anchor_policy=policy,
           seed=42,
           engine='nx'
       )
       # Show anchor index for first few base nodes
       anchors = [gog_p.anchor_vertex(i) for i in range(3)]
       print(f"  Policy '{policy}': anchors = {anchors}")

   # =========================================================
   # Example 3: Mixed graph types
   # =========================================================
   print("\nExample 3: Mixed Graph Types")
   print("-" * 40)

   # 2D lattice base with Erdos-Renyi fibers
   gog_mixed = GraphOfGraphs(
       base_graph_type='Lattice2D',
       base_params={'side1': 4, 'geo': 'sqr'},
       fiber_graph_type='ErdosRenyi',
       fiber_params={'n': 20, 'p': 0.15, 'extract_giant_component': False},
       anchor_policy='first',
       pflip=0.2,
       seed=42,
       engine='nx'
   )

   print(f"Lattice base + ER fibers: {gog_mixed.N} total nodes")
   print(f"Structure: {gog_mixed.N_base} base + {gog_mixed.N_base} x {gog_mixed.N_fiber} fibers")

   # =========================================================
   # Example 4: Accessing structural information
   # =========================================================
   print("\nExample 4: Structural Access")
   print("-" * 40)

   # Get vertex indices for different parts
   base_verts = gog.base_vertex_indices()
   fiber_0_verts = gog.fiber_vertex_indices(0)  # Fiber attached to base node 0
   fiber_5_verts = gog.fiber_vertex_indices(5)  # Fiber attached to base node 5

   print(f"Base vertices: {base_verts[:5]}... ({len(base_verts)} total)")
   print(f"Fiber 0 vertices: {fiber_0_verts[:5]}... ({len(fiber_0_verts)} total)")
   print(f"Fiber 5 vertices: {fiber_5_verts[:5]}... ({len(fiber_5_verts)} total)")

   # Determine which layer a vertex belongs to
   for v in [0, gog.N_base - 1, gog.N_base, gog.N - 1]:
       layer = gog.vertex_layer(v)
       print(f"  Vertex {v}: {layer}")

   # =========================================================
   # Example 5: Spectral analysis
   # =========================================================
   print("\nExample 5: Spectral Analysis")
   print("-" * 40)

   # Check if separated spectrum computation is valid
   can_separate = gog.can_use_separated_spectrum()
   print(f"Can use separated spectrum: {can_separate}")

   if can_separate:
       # Efficient separated computation
       # Returns N_base * N_fiber eigenvalues (Dirac formula)
       separated_eigvals = gog.compute_separated_spectrum()
       print(f"Separated eigenvalues: {len(separated_eigvals)}")
       print(f"Eigenvalue range: [{separated_eigvals.min():.4f}, {separated_eigvals.max():.4f}]")

   # Full spectrum for comparison
   gog.compute_laplacian_spectrum()
   print(f"Full eigenvalues: {len(gog.eigv)}")
   print(f"Full range: [{gog.eigv.min():.4f}, {gog.eigv.max():.4f}]")

   # =========================================================
   # Example 6: Graph-tool backend
   # =========================================================
   print("\nExample 6: Multi-Backend Support")
   print("-" * 40)

   # Same structure with graph-tool backend
   try:
       gog_gt = GraphOfGraphs(
           base_graph_type='Lattice2D',
           base_params={'side1': 5},
           fiber_graph_type='Lattice2D',
           fiber_params={'side1': 4},
           anchor_policy='first',
           seed=42,
           engine='gt'
       )
       print(f"graph-tool backend: {gog_gt.N} nodes")
       print(f"NX and GT produce same structure: N={gog.N} == {gog_gt.N}")
   except ImportError:
       print("graph-tool not available, skipping GT backend example")

   # =========================================================
   # Example 7: Visualization
   # =========================================================
   print("\nExample 7: Visualization")
   print("-" * 40)

   # Create a small graph for visualization
   gog_small = GraphOfGraphs(
       base_graph_type='Lattice2D',
       base_params={'side1': 3},
       fiber_graph_type='Lattice2D',
       fiber_params={'side1': 2},
       anchor_policy='first',
       seed=42,
       engine='nx'
   )

   # Compute spectrum
   gog_small.compute_laplacian_spectrum()

   fig, axes = plt.subplots(1, 2, figsize=(12, 5))

   # Eigenvalue histogram
   axes[0].hist(gog_small.eigv, bins=30, edgecolor='black', alpha=0.7)
   axes[0].axvline(x=0, color='r', linestyle='--', label='Zero')
   axes[0].set_xlabel('Eigenvalue')
   axes[0].set_ylabel('Count')
   axes[0].set_title(f'GraphOfGraphs Spectrum ({gog_small.N} nodes)')
   axes[0].legend()

   # Sorted eigenvalues
   axes[1].plot(np.sort(gog_small.eigv), 'b-', linewidth=1)
   axes[1].axhline(y=0, color='r', linestyle='--', alpha=0.5)
   axes[1].set_xlabel('Index')
   axes[1].set_ylabel('Eigenvalue')
   axes[1].set_title('Sorted Eigenvalues')

   plt.tight_layout()
   plt.savefig('graph_of_graphs_spectrum.png', dpi=150)
   plt.show()

   print(f"\nMetadata: {gog_small.gog_structure.keys()}")

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
   - :doc:`advanced_graphs` - Hierarchical and multi-engine graph structures
   - :doc:`plotting` - Visualization techniques
