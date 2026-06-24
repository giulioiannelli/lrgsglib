Laplacian Renormalization Group
================================

This section describes the theoretical framework of the Laplacian Renormalization
Group (LRG) for signed graphs.

Introduction
------------

The Laplacian Renormalization Group provides a framework for analyzing how
spectral properties of networks change under coarse-graining transformations.
This approach connects network structure to dynamical processes through
spectral analysis.

Key Concepts
~~~~~~~~~~~~

- **Spectral coarse-graining** - Reducing graph complexity while preserving key spectral properties
- **Renormalization flow** - Evolution of spectral characteristics under scale transformation
- **Fixed points** - Network structures that are invariant under renormalization
- **Universality classes** - Groups of networks with similar renormalization behavior

The Signed Laplacian
--------------------

For a signed graph G with N nodes, the **signed Laplacian** is defined as:

.. math::

   L = D - A

where:

- :math:`A` is the signed adjacency matrix with :math:`A_{ij} = \pm 1` for connected nodes
- :math:`D` is the degree matrix with :math:`D_{ii} = \sum_j |A_{ij}|`

Unlike the standard Laplacian, the signed Laplacian can have **negative eigenvalues**
when the graph contains frustration.

Properties
~~~~~~~~~~

1. **Symmetry**: :math:`L = L^T` (for undirected graphs)
2. **Row sum**: :math:`\sum_j L_{ij} = D_{ii} - \sum_j A_{ij}`
3. **Eigenvalue bounds**: For unfrustrated graphs, :math:`0 \leq \lambda_i \leq 2D_{max}`

Spectral Decomposition
----------------------

The eigendecomposition of the signed Laplacian:

.. math::

   L = V \Lambda V^T

where:

- :math:`\Lambda = \text{diag}(\lambda_1, \lambda_2, \ldots, \lambda_N)` are eigenvalues
- :math:`V = [v_1 | v_2 | \ldots | v_N]` are orthonormal eigenvectors

Eigenvalue Ordering
~~~~~~~~~~~~~~~~~~~

Eigenvalues are typically ordered by magnitude:

.. math::

   \lambda_1 \leq \lambda_2 \leq \cdots \leq \lambda_N

For signed graphs, :math:`\lambda_1` may be negative (frustrated) or zero (balanced).

Spectral Gap
~~~~~~~~~~~~

The **spectral gap** :math:`\Delta = \lambda_2 - \lambda_1` characterizes:

- Mixing time of random walks
- Relaxation time of diffusion processes
- Community structure strength

Density Matrix Formalism
------------------------

The spectral density matrix at time scale :math:`\tau`:

.. math::

   \rho(\tau) = \frac{e^{-\tau L}}{\text{Tr}[e^{-\tau L}]}

This density matrix interpolates between:

- :math:`\tau \to 0`: Uniform distribution (identity/N)
- :math:`\tau \to \infty`: Ground state projection

Partition Function
~~~~~~~~~~~~~~~~~~

The partition function:

.. math::

   Z(\tau) = \text{Tr}[e^{-\tau L}] = \sum_{i=1}^{N} e^{-\tau \lambda_i}

encodes all spectral information of the graph.

Information Theory from Spectra
-------------------------------

Von Neumann Entropy
~~~~~~~~~~~~~~~~~~~

The **spectral entropy** at scale :math:`\tau`:

.. math::

   S(\tau) = -\text{Tr}[\rho(\tau) \log \rho(\tau)]

Computed from eigenvalues:

.. math::

   S(\tau) = -\sum_{i=1}^{N} p_i(\tau) \log p_i(\tau)

where :math:`p_i(\tau) = e^{-\tau \lambda_i} / Z(\tau)`.

Specific Heat
~~~~~~~~~~~~~

The **spectral specific heat**:

.. math::

   C(\tau) = \frac{dS}{d\log\tau}

measures the rate of information change across scales.

Spectral Dimension
~~~~~~~~~~~~~~~~~~

The **spectral dimension** :math:`d_s` characterizes how entropy scales:

.. math::

   S(\tau) \sim \frac{d_s}{2} \log \tau \quad \text{as } \tau \to \infty

For regular lattices, :math:`d_s` equals the embedding dimension.
Hierarchical networks often have non-integer :math:`d_s`.

Coarse-Graining Procedure
-------------------------

Basic Algorithm
~~~~~~~~~~~~~~~

1. **Compute spectrum**: Eigendecompose L
2. **Select scale**: Choose number of modes k to retain
3. **Project**: Map nodes to k-dimensional spectral embedding
4. **Cluster**: Group nodes by spectral proximity
5. **Contract**: Create coarse-grained graph
6. **Iterate**: Repeat until desired resolution

Spectral Embedding
~~~~~~~~~~~~~~~~~~

Map node i to k-dimensional vector:

.. math::

   \phi_i = (\sqrt{\lambda_1} v_{1,i}, \sqrt{\lambda_2} v_{2,i}, \ldots, \sqrt{\lambda_k} v_{k,i})

Nodes close in this embedding have similar spectral properties.

Graph Contraction
~~~~~~~~~~~~~~~~~

After clustering, the coarse-grained adjacency:

.. math::

   \tilde{A}_{IJ} = \sum_{i \in I} \sum_{j \in J} A_{ij}

where I, J are clusters.

Renormalization Flow
--------------------

As the graph is coarse-grained through iterations, track:

1. **Spectral gap**: :math:`\Delta^{(n)}` at iteration n
2. **Entropy**: :math:`S^{(n)}(\tau)` profile
3. **Spectral dimension**: :math:`d_s^{(n)}`

Fixed Points
~~~~~~~~~~~~

A graph is at a **fixed point** if its spectral properties are invariant
under coarse-graining:

- **Regular lattices**: Fixed points with :math:`d_s = d`
- **Fractals**: Fixed points with fractional :math:`d_s`
- **Hierarchical**: Self-similar spectral structure

Universality
~~~~~~~~~~~~

Different networks may flow to the same fixed point, indicating shared
large-scale behavior despite microscopic differences.

Applications
------------

Anomalous Diffusion
~~~~~~~~~~~~~~~~~~~

The return probability of a random walker:

.. math::

   P(t) = \frac{1}{N} \text{Tr}[e^{-tL}] = \frac{Z(t)}{N}

For networks with spectral dimension :math:`d_s`:

.. math::

   P(t) \sim t^{-d_s/2}

Anderson Localization
~~~~~~~~~~~~~~~~~~~~~

On frustrated networks, eigenvectors may localize:

.. math::

   \text{IPR}_i = \sum_j v_{i,j}^4

High IPR indicates localization; spectral properties reflect this transition.

Community Detection
~~~~~~~~~~~~~~~~~~~

Spectral methods for community detection use the coarse-graining ideas:

1. Compute low-lying eigenvectors
2. Embed nodes in spectral space
3. Cluster in embedding space

Implementation in lrgsglib
--------------------------

.. code-block:: python

   from lrgsglib.graphs import Lattice2D
   from lrgsglib.utils.lrg.spectral import get_graph_lspectrum
   from lrgsglib.utils.lrg.infocomm import compute_entropy_observables_from_eigenvalues

   # Create frustrated lattice
   lattice = Lattice2D(side=32, pflip=0.3, seed=42)
   lattice.flip_random_fract_edges()

   # Compute full spectrum
   L, eigenvalues = get_graph_lspectrum(lattice.gr['G'], library='scipy')

   # Compute entropy observables
   entropy, heat, variance, time_grid = compute_entropy_observables_from_eigenvalues(
       eigenvalues=eigenvalues,
       num_nodes=lattice.N
   )

   # Spectral dimension from specific heat
   d_s = 2 * heat[-1]  # Asymptotic value

See Also
--------

- :doc:`signed_graphs` - Signed graph theory
- :doc:`stat_phys` - Statistical physics models
- :doc:`/user_guide/spectral` - Practical spectral analysis guide
