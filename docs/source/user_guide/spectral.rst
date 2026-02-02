Spectral Analysis
=================

This guide covers spectral analysis of signed graphs, the core computational
framework of lrgsglib. Spectral properties—eigenvalues and eigenvectors of the
signed Laplacian—reveal fundamental structural and dynamical properties of
networks.

.. contents:: Contents
   :local:
   :depth: 2

Overview
--------

The **signed Laplacian** is defined as:

.. math::

   L_s = D - A_s

where :math:`D` is the degree matrix and :math:`A_s` is the signed adjacency matrix
(with entries :math:`\pm 1`). Unlike unsigned graphs, signed graphs can have
**negative eigenvalues**, which are signatures of frustration—the inability to
satisfy all edge constraints simultaneously.

Key concepts:

- **Spectrum**: The set of eigenvalues :math:`\{\lambda_1, \lambda_2, ..., \lambda_N\}`
- **Spectral gap**: :math:`\lambda_2 - \lambda_1`, related to mixing time
- **Frustration**: Indicated by negative eigenvalues
- **Full spectrum**: Required for entropy calculations (use ``howmany=0``)

Computing the Spectrum
----------------------

lrgsglib provides two main approaches to compute the Laplacian spectrum:

1. **Utility functions** in ``utils.lrg.spectral`` (for quick computations)
2. **SignedGraph methods** (for class-based workflows with caching)

Using Utility Functions
~~~~~~~~~~~~~~~~~~~~~~~

The simplest way to compute a spectrum:

.. code-block:: python

   import networkx as nx
   from lrgsglib.utils.lrg.spectral import get_graph_lspectrum

   # Create a simple graph
   G = nx.path_graph(10)

   # Compute Laplacian and spectrum
   L, eigenvalues = get_graph_lspectrum(G, library='numpy')

   print(f"Laplacian shape: {L.shape}")
   print(f"Eigenvalues: {eigenvalues}")
   print(f"Smallest eigenvalue: {eigenvalues.min():.6f}")

**Available backends:**

- ``'numpy'`` (default): Dense computation, full spectrum
- ``'scipy'`` or ``'sp'``: SciPy's eigenvalue solver
- ``'networkx'`` or ``'nx'``: NetworkX's built-in function
- ``'cupy'`` or ``'cp'``: GPU-accelerated (requires CuPy)

Signed vs Unsigned Laplacian
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, ``get_graph_lspectrum()`` uses NetworkX's unsigned Laplacian, which
ignores edge signs (weights are converted to absolute values). For graphs with
negative edge weights (antiferromagnetic/repulsive interactions), use the
``signed=True`` parameter:

.. code-block:: python

   import networkx as nx
   from lrgsglib.utils.lrg.spectral import get_graph_lspectrum

   # Create a graph with negative edges
   G = nx.Graph()
   G.add_edge(0, 1, weight=-1)  # Antiferromagnetic edge
   G.add_edge(1, 2, weight=1)
   G.add_edge(2, 0, weight=1)

   # Unsigned Laplacian (ignores signs)
   L_unsigned, eigv_unsigned = get_graph_lspectrum(G, signed=False)

   # Signed Laplacian (preserves signs)
   L_signed, eigv_signed = get_graph_lspectrum(G, signed=True)

   print(f"Unsigned eigenvalues: {sorted(eigv_unsigned.real)}")
   print(f"Signed eigenvalues: {sorted(eigv_signed.real)}")

The signed Laplacian is computed as :math:`L_s = D - A_s` where:

- :math:`D` is the degree matrix using **absolute** edge weights
- :math:`A_s` preserves the original edge signs

This convention ensures that the signed Laplacian correctly captures frustration
in antiferromagnetic systems while maintaining proper graph connectivity.

.. note::

   The ``signed_laplacian_matrix()`` function from ``nx_patches.funcs.spectral``
   is the canonical implementation for signed Laplacian matrix construction.
   ``get_graph_lspectrum(..., signed=True)`` uses this implementation internally.

Using SignedGraph Methods
~~~~~~~~~~~~~~~~~~~~~~~~~

For class-based workflows with caching and more control:

.. code-block:: python

   from lrgsglib.nx_patches.Lattice2D import Lattice2D

   # Create a 2D lattice with signed edges
   lattice = Lattice2D(side1=16, geo='sqr', pflip=0.3, seed=42)
   lattice.flip_random_fract_edges()

   # Compute eigenvalues only (faster, less memory)
   lattice.compute_laplacian_spectrum()
   print(f"Eigenvalues computed: {len(lattice.eigv)}")
   print(f"Spectral gap: {lattice.eigv[1] - lattice.eigv[0]:.6f}")

   # Compute full eigendecomposition (eigenvalues + eigenvectors)
   lattice.compute_laplacian_spectrum_weigV()
   print(f"Eigenvectors shape: {lattice.eigV.shape}")

**Key methods:**

- ``compute_laplacian_spectrum()``: Eigenvalues only (fast, memory-efficient)
- ``compute_laplacian_spectrum_weigV()``: Full eigendecomposition
- ``compute_k_eigvV(k=10)``: Compute k smallest eigenvectors

Backend Selection
~~~~~~~~~~~~~~~~~

Choose the backend based on your needs:

.. code-block:: python

   from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph

   # Create a cascade graph
   mc = MultiplicativeCascadeGraph(p1=0.8, p2=0.6, fraction=0.4, seed=42)

   # Force NumPy backend (dense, full spectrum)
   mc.compute_laplacian_spectrum(backend='numpy', keep_sparse=False)

   # Force sparse solver (for large graphs, gets N-2 eigenvalues)
   mc.compute_laplacian_spectrum(backend='scipy', keep_sparse=True)

   # Auto-decide (default): based on graph size and sparsity
   mc.compute_laplacian_spectrum()

.. note::

   Sparse solvers (ARPACK) compute only N-2 eigenvalues. For entropy calculations
   requiring the full spectrum, use ``keep_sparse=False`` or the dense backend.

Entropy and Information Theory
------------------------------

The **von Neumann entropy** measures information content of the graph structure.
For a diffusion process at time :math:`\tau`:

.. math::

   S(\tau) = -\text{Tr}[\rho(\tau) \log \rho(\tau)]

where :math:`\rho(\tau) = \exp(-\tau L) / \text{Tr}[\exp(-\tau L)]` is the
density matrix.

Computing Entropy from Eigenvalues
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The standard approach computes entropy from the full spectrum:

.. code-block:: python

   from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph
   from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
   from lrgsglib.utils.lrg.infocomm import compute_entropy_observables_from_eigenvalues

   # Create graph
   mc = MultiplicativeCascadeGraph(p1=0.8, p2=0.6, fraction=0.4, seed=42)

   # Compute full spectrum
   L, eigenvalues = get_graph_lspectrum(mc.gr['G'], library='numpy')

   # Compute entropy observables
   entropy, specific_heat, variance, time_grid = compute_entropy_observables_from_eigenvalues(
       eigenvalues=eigenvalues,
       num_nodes=mc.N,
       steps=200,          # Number of time points
       t1=-2, t2=5,        # log10 time range
       entropy_norm='complement',  # Returns 1 - S/log(N)
   )

   print(f"Time range: [{time_grid[0]:.2e}, {time_grid[-1]:.2e}]")
   print(f"Entropy at tau=1: {entropy[100]:.4f}")

**Output observables:**

- ``entropy``: Normalized von Neumann entropy profile
- ``specific_heat``: Derivative :math:`C(\tau) = dS/d(\log\tau)` (spectral dimension proxy)
- ``variance``: Variance of eigenvalue distribution
- ``time_grid``: Logarithmic time points

Matrix-Free Entropy (Large Graphs)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For very large graphs where eigendecomposition is expensive, use matrix-free methods:

.. code-block:: python

   from lrgsglib.utils.lrg.infocomm import (
       compute_entropy_observables_slq,        # Stochastic Lanczos Quadrature
       compute_entropy_observables_expm_multiply,  # Matrix exponential
   )

   # For large sparse Laplacians
   mc.upd_graph_matrices()  # Ensure slp (signed Laplacian) is computed

   # Method 1: Stochastic Lanczos Quadrature (SLQ)
   entropy, heat, var, tau = compute_entropy_observables_slq(
       laplacian=mc.slp,
       num_nodes=mc.N,
       steps=100,
       lanczos_steps=50,
       num_samples=30,
       seed=42,
   )

   # Method 2: expm_multiply (matrix exponential action)
   entropy, heat, var, tau = compute_entropy_observables_expm_multiply(
       L=mc.slp,
       num_nodes=mc.N,
       steps=100,
       num_samples=30,
       seed=42,
   )

.. warning::

   Matrix-free methods are approximations. Use eigenvalue-based computation for
   publication-quality results on small-to-medium graphs (N < 50,000).

Spectral Dimension
~~~~~~~~~~~~~~~~~~

The **spectral dimension** characterizes how diffusion scales on the graph:

.. math::

   d_s \approx 2 \langle C(\tau \to \infty) \rangle

Extract it from the specific heat tail:

.. code-block:: python

   import numpy as np

   # Use tail 10% of specific heat profile
   tail_fraction = 0.1
   tail_len = int(tail_fraction * len(specific_heat))
   ds_estimate = 2.0 * np.mean(specific_heat[-tail_len:])
   print(f"Spectral dimension estimate: {ds_estimate:.3f}")

Renyi and Tsallis Entropies
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generalized entropies provide additional structural information:

.. code-block:: python

   from lrgsglib.utils.lrg.infocomm import compute_renyi_observables_from_eigenvalues

   # Compute Renyi entropy of order q=2
   results = compute_renyi_observables_from_eigenvalues(
       eigenvalues=eigenvalues,
       num_nodes=mc.N,
       q=2.0,
       steps=200,
   )

   print(f"Renyi entropy at tau=1: {results['renyi_entropy'][100]:.4f}")
   print(f"Tsallis entropy at tau=1: {results['tsallis_entropy'][100]:.4f}")
   print(f"Spectral dimension (ds): {results['ds_estimate']:.3f}")

Eigenvector Analysis
--------------------

Eigenvectors reveal community structure and localization patterns.

Computing Eigenvectors
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from lrgsglib.nx_patches.Lattice2D import Lattice2D

   # Create frustrated lattice
   lattice = Lattice2D(side1=32, geo='sqr', pflip=0.2, seed=42)
   lattice.flip_random_fract_edges()

   # Compute full eigendecomposition
   lattice.compute_laplacian_spectrum_weigV(
       transpose=True,      # eigV[i] gives i-th eigenvector
       flip_to_pos=True,    # Normalize sign for consistency
   )

   # Access eigenvectors
   fiedler_vector = lattice.eigV[1]  # Second eigenvector (Fiedler)
   print(f"Fiedler vector shape: {fiedler_vector.shape}")
   print(f"Fiedler eigenvalue: {lattice.eigv[1]:.6f}")

   # Compute k smallest eigenvectors only (more efficient)
   lattice.compute_k_eigvV(k=5)

Eigenvector-Based Clustering
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Fiedler vector (second eigenvector) partitions the graph:

.. code-block:: python

   import numpy as np

   # Get Fiedler vector
   fiedler = lattice.eigV[1]

   # Binary partition based on sign
   cluster1 = np.where(fiedler > 0)[0]
   cluster2 = np.where(fiedler <= 0)[0]

   print(f"Cluster 1 size: {len(cluster1)}")
   print(f"Cluster 2 size: {len(cluster2)}")

For more sophisticated clustering, use the ``get_eigV_cluster_sizes()`` method
from SignedGraph's partitioning module.

Spectral Embeddings
~~~~~~~~~~~~~~~~~~~

Project nodes into a low-dimensional spectral space:

.. code-block:: python

   # Use first k eigenvectors as coordinates
   k = 3
   embedding = lattice.eigV[:k].T  # Shape: (N, k)

   print(f"Spectral embedding shape: {embedding.shape}")
   print(f"First node coordinates: {embedding[0]}")

Frustration Analysis
--------------------

Negative eigenvalues indicate frustration in signed graphs.

Detecting Frustration
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from lrgsglib.nx_patches.ErdosRenyi import ErdosRenyi

   # Create signed Erdos-Renyi graph
   er = ErdosRenyi(nnodes=100, prob=0.1, pflip=0.3, seed=42)
   er.flip_random_fract_edges()

   # Compute spectrum
   er.compute_laplacian_spectrum()

   # Count negative eigenvalues
   negative_eigvals = er.eigv[er.eigv < 0]
   print(f"Number of negative eigenvalues: {len(negative_eigvals)}")
   print(f"Most negative eigenvalue: {er.eigv.min():.4f}")
   print(f"Frustration indicator: {len(negative_eigvals) / er.N:.2%}")

Frustration and Balance
~~~~~~~~~~~~~~~~~~~~~~~

A signed graph is **balanced** if it can be partitioned into two groups where:

- All positive edges are within groups
- All negative edges are between groups

Frustration prevents perfect balance:

.. code-block:: python

   # Spectral gap indicates balance
   # Small gap near zero → high frustration
   spectral_gap = er.eigv[1] - er.eigv[0]
   print(f"Spectral gap: {spectral_gap:.4f}")

   # Large negative eigenvalue magnitude → more frustration
   max_frustration = abs(er.eigv.min())
   print(f"Max frustration magnitude: {max_frustration:.4f}")

Working with Different Graph Types
----------------------------------

MultiplicativeCascade (Hierarchical)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph

   mc = MultiplicativeCascadeGraph(
       p1=0.8, p2=0.6, p3=0.6, p4=0.8,
       fraction=0.4,
       iterations=6,
       seed=42,
   )

   # Access the NetworkX graph for spectral functions
   G = mc.gr['G']  # Integer-labeled nodes
   # or
   H = mc.gr['H']  # Coordinate-labeled nodes (for lattices)

   # Class method
   mc.compute_laplacian_spectrum()

Lattice2D/3D (Regular Grids)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from lrgsglib.nx_patches.Lattice2D import Lattice2D
   from lrgsglib.nx_patches.Lattice3D import Lattice3D

   # 2D square lattice
   l2d = Lattice2D(side1=32, geo='sqr', pflip=0.2, seed=42)
   l2d.flip_random_fract_edges()
   l2d.compute_laplacian_spectrum()

   # 3D cubic lattice
   l3d = Lattice3D(side1=10, geo='cub', pflip=0.1, seed=42)
   l3d.flip_random_fract_edges()
   l3d.compute_laplacian_spectrum()

Best Practices
--------------

1. **Always set seeds** for reproducibility:

   .. code-block:: python

      graph = Lattice2D(side1=32, pflip=0.3, seed=42)

2. **Use eigenvalues-only** when eigenvectors aren't needed:

   .. code-block:: python

      # Fast (eigenvalues only)
      graph.compute_laplacian_spectrum()

      # Slower (full decomposition)
      graph.compute_laplacian_spectrum_weigV()

3. **Choose backends wisely**:

   - Small graphs (N < 5000): ``'numpy'`` (dense, exact)
   - Large sparse graphs: ``'scipy'`` with ``keep_sparse=True``
   - GPU available: ``'cupy'`` for maximum performance

4. **Full spectrum for entropy**: Use ``keep_sparse=False`` or ``howmany=0``

5. **Call** ``flip_random_fract_edges()`` **after creating lattices with** ``pflip``:

   .. code-block:: python

      lattice = Lattice2D(side1=32, pflip=0.3)
      lattice.flip_random_fract_edges()  # Actually applies the flips!

API Reference
-------------

.. seealso::

   - :mod:`lrgsglib.utils.lrg.spectral` - Utility functions
   - :mod:`lrgsglib.utils.lrg.infocomm` - Entropy computations
   - :class:`lrgsglib.nx_patches.SignedGraph.SignedGraph` - Base class methods
