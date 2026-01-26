Advanced Graph Types
====================

This section covers advanced graph structures with specialized hierarchical
and fractal properties, including multiplicative cascades, Vicsek graphs,
and Dirac lattice structures.

.. contents:: On This Page
   :local:
   :depth: 2

Overview
--------

Beyond standard lattices and random graphs, lrgsglib provides specialized
graph generators for studying:

- **Hierarchical structures** via multiplicative cascades
- **Fractal graphs** via Kronecker products (Vicsek)
- **Product structures** via Dirac combs and brushes

These graph types inherit from :class:`MultispectralGraph`, providing
efficient spectral computation methods.

Multiplicative Cascade Graphs
-----------------------------

Multiplicative cascades generate hierarchical networks with tunable
multifractal properties.

Basic Usage
^^^^^^^^^^^

.. code-block:: python

   from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph

   # Create cascade with seed probabilities
   mc = MultiplicativeCascadeGraph(
       p1=0.8, p2=0.6, p3=0.6, p4=0.4,
       fraction=0.5,      # Sample 50% of nodes
       iterations=8,      # Cascade depth
       variant="exp_clocks"  # Fast algorithm (default)
   )

   print(f"Nodes: {mc.N}, Edges: {mc.Ne}")
   print(f"Probability matrix shape: {mc.probability_matrix.shape}")

Parameters
^^^^^^^^^^

The cascade is controlled by four seed probabilities that define the
hierarchical edge probabilities:

+-------------+----------------------------------------------+
| Parameter   | Description                                  |
+=============+==============================================+
| ``p1``      | Top-left quadrant probability (default: 0.9) |
+-------------+----------------------------------------------+
| ``p2``      | Top-right quadrant probability (default: 0.7)|
+-------------+----------------------------------------------+
| ``p3``      | Bottom-left quadrant probability (default: 0.7)|
+-------------+----------------------------------------------+
| ``p4``      | Bottom-right quadrant probability (default: 0.5)|
+-------------+----------------------------------------------+
| ``fraction``| Fraction of nodes to sample (default: 1.0)   |
+-------------+----------------------------------------------+
| ``iterations``| Number of cascade iterations (default: 10) |
+-------------+----------------------------------------------+
| ``variant`` | Algorithm: "exp_clocks" (fast) or "standard" |
+-------------+----------------------------------------------+

Algorithm Variants
^^^^^^^^^^^^^^^^^^

Two generation algorithms are available:

**exp_clocks (default):**
  Uses exponential clocks for 2-16x faster generation. Recommended for
  most use cases.

**standard:**
  Original algorithm, useful for validation or when exact reproducibility
  with older code is needed.

.. code-block:: python

   # Compare variants
   mc_fast = MultiplicativeCascadeGraph(
       p1=0.8, p2=0.6, p3=0.6, p4=0.4,
       fraction=0.3, iterations=10,
       variant="exp_clocks"
   )

   mc_standard = MultiplicativeCascadeGraph(
       p1=0.8, p2=0.6, p3=0.6, p4=0.4,
       fraction=0.3, iterations=10,
       variant="standard"
   )

Spectral Analysis
^^^^^^^^^^^^^^^^^

Multiplicative cascades exhibit characteristic spectral properties:

.. code-block:: python

   from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph
   from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
   from lrgsglib.utils.lrg.infocomm import compute_entropy_observables_from_eigenvalues

   # Create cascade
   mc = MultiplicativeCascadeGraph(
       p1=0.8, p2=0.6, p3=0.6, p4=0.4,
       fraction=0.4, iterations=8
   )

   # Compute spectrum
   L, eigenvalues = get_graph_lspectrum(mc.gr['G'], library='numpy')

   # Compute entropy observables
   entropy, heat, var, tau = compute_entropy_observables_from_eigenvalues(
       eigenvalues=eigenvalues,
       num_nodes=mc.N
   )

   print(f"Spectral dimension estimate: {2 * heat[-1]:.3f}")


Vicsek Graphs
-------------

Vicsek graphs (Palla-Lovász-Vicsek) are fractal networks generated through
iterated Kronecker products, exhibiting exact self-similarity.

Basic Usage
^^^^^^^^^^^

.. code-block:: python

   from lrgsglib.nx_patches.Vicsek import VicsekGraph

   # Create Vicsek graph
   vg = VicsekGraph(
       N=1000,        # Target number of nodes
       k=4,           # Kronecker iterations
       m=3,           # Initial measure dimension
       symmetric=True # Symmetrize initial measure
   )

   print(f"Nodes: {vg.N}, Edges: {vg.Ne}")

Parameters
^^^^^^^^^^

+-------------+----------------------------------------------+
| Parameter   | Description                                  |
+=============+==============================================+
| ``N``       | Target number of nodes                       |
+-------------+----------------------------------------------+
| ``k``       | Number of Kronecker product iterations       |
+-------------+----------------------------------------------+
| ``pij``     | Custom initial measure matrix (optional)     |
+-------------+----------------------------------------------+
| ``m``       | Dimension of auto-generated measure matrix   |
+-------------+----------------------------------------------+
| ``symmetric``| Symmetrize the initial measure             |
+-------------+----------------------------------------------+

Custom Initial Measure
^^^^^^^^^^^^^^^^^^^^^^

You can provide a custom initial measure matrix:

.. code-block:: python

   import numpy as np
   from lrgsglib.nx_patches.Vicsek import VicsekGraph

   # Define custom initial measure
   pij = np.array([
       [0.8, 0.2, 0.1],
       [0.2, 0.9, 0.3],
       [0.1, 0.3, 0.7]
   ])

   vg = VicsekGraph(N=500, k=3, pij=pij)
   print(f"Generated {vg.N} nodes with custom measure")

Mathematical Background
^^^^^^^^^^^^^^^^^^^^^^^

The Vicsek graph is constructed by taking iterated Kronecker products
of an initial measure matrix. After k iterations, the probability of
an edge between nodes i and j is given by:

.. math::

   P_{ij}^{(k)} = \prod_{l=1}^{k} p_{i_l j_l}

where :math:`(i_1, ..., i_k)` and :math:`(j_1, ..., j_k)` are the
k-ary representations of i and j respectively.


Dirac Lattice Structures
------------------------

Dirac lattices are product structures consisting of a base graph with
fiber graphs attached to each base node. This structure enables
efficient spectral computation by separating the computation into
base and fiber components.

Dirac Comb (1D Base + 1D Fibers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A Dirac comb has a 1D chain as the base graph with 1D chains as fibers:

.. code-block:: python

   from lrgsglib.nx_patches.DiracLattice import DiracCombGraph

   # Create Dirac comb
   comb = DiracCombGraph(
       base_nodes=20,     # 1D base with 20 nodes
       fiber_nodes=10,    # Each base node has 10-node fiber
       periodic=True      # Periodic boundary conditions
   )

   print(f"Total nodes: {comb.N}")  # 20 * (1 + 10) = 220
   print(f"Base dimension: {comb.base_nodes}")
   print(f"Fiber dimension: {comb.fiber_nodes}")

Dirac Brush (2D Base + 1D Fibers)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A Dirac brush has a 2D grid as the base graph with 1D chains as fibers:

.. code-block:: python

   from lrgsglib.nx_patches.DiracLattice import DiracBrushGraph

   # Create Dirac brush
   brush = DiracBrushGraph(
       base_x=10,         # 10x10 2D base
       base_y=10,
       fiber_nodes=5,     # 5-node fiber at each base point
       periodic=True
   )

   print(f"Total nodes: {brush.N}")  # 100 * (1 + 5) = 600
   print(f"Base shape: {brush.base_x}x{brush.base_y}")

Efficient Spectral Computation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The key advantage of Dirac structures is efficient spectral computation.
Instead of diagonalizing the full N×N Laplacian, we can exploit the
product structure:

**Complexity comparison:**

- Full diagonalization: O((N_base × N_fiber)³)
- Separated computation: O(N_base³ + N_base × N_fiber³)

For a brush with 100 base nodes and 10 fiber nodes (1000 total):

- Full: O(10⁹) operations
- Separated: O(10⁶ + 10⁵) = O(10⁶) operations

.. code-block:: python

   from lrgsglib.nx_patches.DiracLattice import DiracBrushGraph
   import numpy as np

   # Create brush
   brush = DiracBrushGraph(
       base_x=20, base_y=20,
       fiber_nodes=8,
       periodic=True
   )

   # Compute spectrum efficiently using separated method
   eigenvalues = brush.compute_dirac_spectrum_separated()

   print(f"Computed {len(eigenvalues)} eigenvalues")
   print(f"Smallest: {eigenvalues[0]:.6f}")
   print(f"Largest: {eigenvalues[-1]:.6f}")

   # Access base and fiber spectra separately
   base_eigenvalues = brush.compute_base_spectrum()
   fiber_lap = brush.compute_fiber_laplacian()

   print(f"Base eigenvalues: {len(base_eigenvalues)}")
   print(f"Fiber Laplacian shape: {fiber_lap.shape}")

Mathematical Background
^^^^^^^^^^^^^^^^^^^^^^^

For a Dirac structure with base eigenvalues :math:`\{\lambda_i^{(b)}\}`
and fiber Laplacian :math:`L_f`, the full spectrum is obtained by solving:

.. math::

   (L_f + \lambda_i^{(b)} e_1 e_1^T) v = \mu v

for each base eigenvalue :math:`\lambda_i^{(b)}`, where :math:`e_1` is
the unit vector corresponding to the attachment point.


Working with Hierarchical Structures
------------------------------------

Common Patterns
^^^^^^^^^^^^^^^

All hierarchical graph types share common patterns:

.. code-block:: python

   from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph
   from lrgsglib.nx_patches.Vicsek import VicsekGraph
   from lrgsglib.nx_patches.DiracLattice import DiracCombGraph

   # All inherit from MultispectralGraph and SignedGraph
   mc = MultiplicativeCascadeGraph(p1=0.8, p2=0.6, iterations=6, fraction=0.5)
   vg = VicsekGraph(N=500, k=3)
   dc = DiracCombGraph(base_nodes=20, fiber_nodes=10)

   # Access underlying NetworkX graph
   for graph_obj in [mc, vg, dc]:
       G = graph_obj.gr['G']  # Integer-labeled graph
       H = graph_obj.gr['H']  # Coordinate-labeled graph (if available)
       print(f"{type(graph_obj).__name__}: N={graph_obj.N}, E={graph_obj.Ne}")

Adding Frustration
^^^^^^^^^^^^^^^^^^

Hierarchical graphs can have signed edges added:

.. code-block:: python

   from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph

   # Create cascade with frustration
   mc = MultiplicativeCascadeGraph(
       p1=0.8, p2=0.6, p3=0.6, p4=0.4,
       fraction=0.4, iterations=7,
       pflip=0.2  # 20% negative edges
   )

   # Apply the edge flips
   mc.flip_random_fract_edges()

   # Count negative edges
   n_neg = sum(1 for u, v in mc.gr['G'].edges()
               if mc.gr['G'][u][v].get('sign', 1) == -1)
   print(f"Negative edges: {n_neg} / {mc.Ne}")


Performance Tips
----------------

Backend Selection
^^^^^^^^^^^^^^^^^

For large hierarchical graphs, choose backends carefully:

.. code-block:: python

   from lrgsglib.utils.lrg.spectral import get_graph_lspectrum

   # For dense graphs (N < 5000)
   L, eigvals = get_graph_lspectrum(mc.gr['G'], library='numpy')

   # For sparse graphs or partial spectrum
   L, eigvals = get_graph_lspectrum(mc.gr['G'], library='scipy')

   # For GPU acceleration (if available)
   L, eigvals = get_graph_lspectrum(mc.gr['G'], library='cupy')

Memory Considerations
^^^^^^^^^^^^^^^^^^^^^

Multiplicative cascades can generate very large graphs. Use the
``fraction`` parameter to control size:

.. code-block:: python

   # Memory-efficient: sample 10% of potential nodes
   mc_small = MultiplicativeCascadeGraph(
       p1=0.8, p2=0.6, p3=0.6, p4=0.4,
       fraction=0.1,  # Small sample
       iterations=12  # Many iterations but sparse
   )

   # Full graph (memory intensive for large iterations)
   mc_full = MultiplicativeCascadeGraph(
       p1=0.8, p2=0.6, p3=0.6, p4=0.4,
       fraction=1.0,
       iterations=8  # Fewer iterations
   )


See Also
--------

- :doc:`graphs` - Basic signed graph operations
- :doc:`lattices` - Regular lattice structures
- :doc:`spectral` - Spectral analysis methods
- :doc:`examples` - Complete workflow examples
