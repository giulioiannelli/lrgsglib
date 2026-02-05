Signed Graph Theory
===================

This section covers the mathematical foundations of signed graphs, balance theory,
and frustration metrics.

Definition
----------

A **signed graph** is a graph G = (V, E) where each edge e ∈ E has an associated
sign σ(e) ∈ {+1, -1}:

- **Positive edges** (+1): Represent agreement, attraction, ferromagnetic coupling
- **Negative edges** (-1): Represent disagreement, repulsion, antiferromagnetic coupling

Signed Adjacency Matrix
~~~~~~~~~~~~~~~~~~~~~~~

The signed adjacency matrix A has entries:

.. math::

   A_{ij} = \begin{cases}
   +1 & \text{if } (i,j) \in E \text{ and } \sigma(i,j) = +1 \\
   -1 & \text{if } (i,j) \in E \text{ and } \sigma(i,j) = -1 \\
   0 & \text{if } (i,j) \notin E
   \end{cases}

Signed Laplacian
~~~~~~~~~~~~~~~~

The **signed Laplacian** matrix:

.. math::

   L = D - A

where D is the absolute degree matrix: :math:`D_{ii} = \sum_j |A_{ij}|`.

Key properties:

- Symmetric for undirected graphs
- Can have negative eigenvalues (unlike standard Laplacian)
- Negative eigenvalues indicate frustration

Balance Theory
--------------

A signed graph is **balanced** if it can be partitioned into two groups such that:

- All positive edges are within groups
- All negative edges are between groups

Characterization Theorem
~~~~~~~~~~~~~~~~~~~~~~~~

A signed graph is balanced if and only if:

1. It contains no cycles with an odd number of negative edges
2. All eigenvalues of the signed Laplacian are non-negative
3. The smallest eigenvalue λ₁ = 0

Balanced Partition
~~~~~~~~~~~~~~~~~~

For a balanced graph, the partition can be found from the eigenvector
corresponding to λ₁ = 0:

.. math::

   v_1 = (s_1, s_2, \ldots, s_N)^T

where :math:`s_i \in \{+1, -1\}` indicates group membership.

Frustration
-----------

**Frustration** occurs when not all constraints can be satisfied simultaneously.
In signed graphs, frustration arises from cycles with an odd number of negative edges.

Frustration Index
~~~~~~~~~~~~~~~~~

The **frustration index** :math:`\phi(G)` is the minimum number of edges that must
be removed (or have their sign changed) to make the graph balanced:

.. math::

   \phi(G) = \min_{\sigma'} |\{e : \sigma(e) \neq \sigma'(e)\}|

where the minimum is over all balanced sign assignments σ'.

Line Index of Balance
~~~~~~~~~~~~~~~~~~~~~

The **line index of balance** L(G) normalizes the frustration index:

.. math::

   L(G) = 1 - \frac{\phi(G)}{|E|}

- L(G) = 1: Perfectly balanced
- L(G) = 0: Maximum frustration (half edges wrong)

Spectral Frustration
~~~~~~~~~~~~~~~~~~~~

The smallest eigenvalue λ₁ indicates frustration:

- λ₁ = 0: Balanced graph
- λ₁ < 0: Frustrated graph
- |λ₁|: Magnitude of frustration

Triangle Frustration
~~~~~~~~~~~~~~~~~~~~

A triangle with vertices i, j, k is **frustrated** if:

.. math::

   A_{ij} \cdot A_{jk} \cdot A_{ki} = -1

The fraction of frustrated triangles is a local frustration measure.

Structural Balance
------------------

**Structural balance** theory (from social network analysis) states:

1. "The friend of my friend is my friend" (+)(+) = (+)
2. "The enemy of my enemy is my friend" (−)(−) = (+)
3. "The friend of my enemy is my enemy" (+)(−) = (−)
4. "The enemy of my friend is my enemy" (−)(+) = (−)

k-Balance
~~~~~~~~~

A generalization allows k ≥ 2 groups:

- All positive edges within groups
- Negative edges between different groups

The graph is **k-balanced** if it can be partitioned into k groups satisfying
these constraints.

Signed Graph Operations
-----------------------

Switching
~~~~~~~~~

**Switching** at node i negates all edges incident to i:

.. math::

   A'_{ij} = \begin{cases}
   -A_{ij} & \text{if } j \text{ is neighbor of } i \\
   A_{ij} & \text{otherwise}
   \end{cases}

Two graphs are **switching equivalent** if one can be obtained from the
other by a sequence of switchings. Switching preserves:

- Balance
- Frustration index
- Eigenvalues of the signed Laplacian

Negation
~~~~~~~~

**Negation** multiplies all edge signs by -1:

.. math::

   A' = -A

This transforms:

- Eigenvalues: λ'ᵢ = -λᵢ + 2D (shifted and negated)
- Frustrated cycles become balanced and vice versa

Signed Random Graphs
--------------------

Signed Erdős-Rényi
~~~~~~~~~~~~~~~~~~

A signed Erdős-Rényi graph G(n, p, q) has:

- n nodes
- Each pair connected with probability p
- Each edge is negative with probability q

Expected properties:

- Number of edges: p·n(n-1)/2
- Fraction negative: q
- Expected frustration depends on q and p

Signed Lattices
~~~~~~~~~~~~~~~

Regular lattices with random sign assignment:

.. code-block:: python

   from lrgsglib.graphs.Lattice2D import Lattice2D

   # Create lattice with 30% negative edges
   lattice = Lattice2D(side=32, pflip=0.3, seed=42)
   lattice.flip_random_fract_edges()

The parameter ``pflip`` controls the fraction of edges flipped to negative.

Spectral Properties
-------------------

Eigenvalue Distribution
~~~~~~~~~~~~~~~~~~~~~~~

For random signed graphs:

- **Balanced** (pflip=0): All eigenvalues ≥ 0, semicircular bulk
- **Frustrated** (pflip>0): Negative eigenvalues appear
- **Maximum frustration** (pflip=0.5): Symmetric distribution around 0

Localization
~~~~~~~~~~~~

On frustrated graphs, eigenvectors may **localize**:

- Low-lying eigenvectors concentrate on frustrated regions
- Inverse Participation Ratio (IPR) measures localization

.. math::

   \text{IPR}_i = \sum_j v_{i,j}^4

- Delocalized: IPR ∼ 1/N
- Localized: IPR ∼ O(1)

Connection to Physics
---------------------

Spin Glasses
~~~~~~~~~~~~

Signed graphs model spin glass systems:

- Nodes = spins (±1)
- Positive edges = ferromagnetic coupling (align)
- Negative edges = antiferromagnetic coupling (anti-align)

Ground state energy minimization becomes balance optimization.

Magnetic Frustration
~~~~~~~~~~~~~~~~~~~~

In magnetic systems, frustration prevents simultaneous satisfaction of
all pairwise interactions, leading to:

- Degenerate ground states
- Slow dynamics
- Complex energy landscapes

Social Networks
~~~~~~~~~~~~~~~

In social network analysis:

- Positive edges = friendship, trust
- Negative edges = enmity, distrust
- Balance = stable social configurations

Implementation in lrgsglib
--------------------------

Creating Signed Graphs
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   from lrgsglib.graphs.SignedGraph import SignedGraph
   from lrgsglib.graphs.ErdosRenyi import signed_erdos_renyi
   import networkx as nx

   # From existing NetworkX graph
   G = nx.erdos_renyi_graph(100, 0.1, seed=42)
   sg = SignedGraph(G, pflip=0.3, seed=42)
   sg.flip_random_fract_edges()

   # Direct construction
   sg = signed_erdos_renyi(n=100, p=0.1, pflip=0.3, seed=42)

Analyzing Frustration
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Compute spectrum
   sg.compute_laplacian_spectrum_weigV()

   # Check for frustration
   lambda_min = sg.eigv[0]
   if lambda_min < -1e-10:
       print(f"Graph is frustrated: λ₁ = {lambda_min:.4f}")
   else:
       print("Graph is balanced")

   # Count frustrated triangles
   n_frustrated = sg.count_frustrated_triangles()

See Also
--------

- :doc:`lrg_theory` - Laplacian renormalization group
- :doc:`stat_phys` - Statistical physics models
- :doc:`/user_guide/graphs` - Practical graph operations
