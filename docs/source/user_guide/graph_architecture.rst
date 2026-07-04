Graph Architecture
==================

This section documents the full graph type system in ``lrgsglib``: the
three-layer architecture (facade, engine, implementation), the inheritance
hierarchy, and the complete catalogue of available graph types with their
parameters and backend support.

.. contents:: On this page
   :local:
   :depth: 2

Three-Layer Architecture
------------------------

Every graph in ``lrgsglib`` follows a **facade / engine / implementation**
pattern:

.. code-block:: text

   Layer 1 — Facade (engine-agnostic)
       lrgsglib.graphs.Lattice2D          ← what users import
       lrgsglib.graphs.ErdosRenyi
       lrgsglib.graphs.TemporalGraph
       ...

   Layer 2 — Engine dispatch (_engine.py)
       GraphEngine.NETWORKX  →  lazy import  →  Lattice2DNX
       GraphEngine.GRAPHTOOL →  lazy import  →  Lattice2DGT

   Layer 3 — Implementation (engine-specific)
       lrgsglib.graphs.nx.Lattice2DNX     ← NX-specific code
       lrgsglib.graphs.gt.Lattice2DGT     ← GT-specific code

**Facades** are classes in ``lrgsglib/graphs/`` whose ``__new__`` method
dispatches to the correct engine implementation. They accept an ``engine=``
parameter (``'nx'``, ``'gt'``, or ``None`` for default) and translate
parameters between engines when needed.

.. code-block:: python

   from lrgsglib.graphs import Lattice2D

   # These return different types but the same logical graph:
   lat_nx = Lattice2D(side1=50, geo='sqr', engine='nx')  # → Lattice2DNX
   lat_gt = Lattice2D(side1=50, geo='sqr', engine='gt')  # → Lattice2DGT

Engine Selection
~~~~~~~~~~~~~~~~

Three ways to choose the backend, in priority order:

1. **Explicit parameter**: ``Lattice2D(..., engine='gt')``
2. **Environment variable**: ``LRGSG_GRAPH_ENGINE=gt``
3. **Global setting**: ``set_default_engine('gt')``

If no engine is specified, NetworkX is the default.

.. code-block:: python

   from lrgsglib.graphs import set_default_engine, get_default_engine

   set_default_engine('gt')   # all subsequent calls use graph-tool
   print(get_default_engine())  # GraphEngine.GRAPHTOOL


Inheritance Hierarchy
---------------------

All concrete graph types ultimately inherit from **SignedGraph** (NX or GT),
which is the structural base providing signed-edge operations, spectral
analysis, entropy computation, clustering, and dynamics integration.

.. code-block:: text

   SignedGraphNX / SignedGraphGT
   ├── Lattice2DNX / Lattice2DGT        (2D lattices, 4+ geometries)
   ├── Lattice3DNX / Lattice3DGT        (3D lattices, 3 geometries)
   ├── LatticeNDNX / LatticeNDGT        (generic cubic, any dimension)
   ├── ErdosRenyiNX / ErdosRenyiGT      (random graphs)
   ├── BarabasiAlbertNX / BarabasiAlbertGT
   ├── WattsStrogatzNX / WattsStrogatzGT
   ├── StochasticBlockModelNX / StochasticBlockModelGT
   ├── kRegularGraphNX / kRegularGraphGT
   ├── ConfigurationModelNX / ...GT
   ├── RandomGeometricNX / ...GT
   ├── ExtendedBarabasiAlbertNX / ...GT
   ├── DualBarabasiAlbertNX / ...GT
   ├── HolmeKimNX / ...GT
   ├── LFRBenchmarkNX / ...GT
   ├── FullyConnectedNX / ...GT
   ├── BipartiteGraphNX / ...GT
   ├── BipartiteFromDegreeSequenceNX / ...GT
   ├── HofieldNNNX / ...GT
   ├── SCSGeneralizedNNNX / ...GT
   ├── MultispectralGraphNX (abstract)
   │   ├── MultiplicativeCascadeNX / ...GT
   │   ├── VicsekNX / ...GT
   │   ├── GraphOfGraphsNX / ...GT
   │   │   ├── DiracCombGraphNX / ...GT
   │   │   └── DiracBrushGraphNX / ...GT
   │   └── HierarchicalModularNX / ...GT
   ├── FractalGraphNX (abstract)
   │   ├── SierpinskiNX / ...GT
   │   └── DGMgraphNX / ...GT
   └── TemporalGraphNX / ...GT          (time-varying, inherits SignedGraphNX)

**Abstract bases** (MultispectralGraphNX, FractalGraphNX, DiracLatticeNX,
CompleteGraphNX, RandomGraphNX) are not directly instantiable — use their
concrete subclasses.


SignedGraph: The Base Class
---------------------------

Every graph in ``lrgsglib`` is a signed graph. The base class
(``SignedGraphNX`` / ``SignedGraphGT``) provides:

Core Attributes
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   * - Attribute
     - Type
     - Description
   * - ``gr``
     - dict
     - Graph representation dictionary. ``gr['G']`` is the primary
       integer-labeled NetworkX graph. Lattices also have ``gr['H']``
       with coordinate-labeled nodes.
   * - ``N``
     - int (property)
     - Number of nodes.
   * - ``Ne``
     - int (property)
     - Number of edges.
   * - ``Ne_n``
     - int (property)
     - Number of negative edges.
   * - ``pflip``
     - float
     - Fraction of edges marked for sign flipping (0.0 to 1.0).
   * - ``seed``
     - int
     - Random seed for reproducibility.
   * - ``on_g``
     - str
     - Active graph representation identifier (default ``'G'``).
   * - ``eset``
     - dict
     - Maps representation → set of all edges.
   * - ``fleset``
     - dict
     - Maps representation → set of negative edges.
   * - ``lfeset``
     - dict
     - Maps representation → set of positive edges.

Matrix Properties (lazy-computed)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   * - Property
     - Type
     - Description
   * - ``adj``
     - scipy.sparse
     - Signed adjacency matrix (weights are ±1).
   * - ``lap``
     - scipy.sparse
     - Standard Laplacian L = D - A.
   * - ``slp``
     - scipy.sparse
     - Signed Laplacian L_s = D_s - A (absolute degrees).
   * - ``degm``
     - scipy.sparse
     - Degree matrix D.
   * - ``sdeg``
     - scipy.sparse
     - Signed degree matrix D_s.
   * - ``eigv``
     - ndarray
     - Eigenvalues of the signed Laplacian.
   * - ``eigV``
     - ndarray
     - Eigenvectors (row-major by default in NX).

Method Categories
~~~~~~~~~~~~~~~~~

**Edge operations:**
  ``flip_random_fract_edges()``, ``flip_sel_edges(links)``,
  ``get_random_links(n)``

**Spectral analysis** (``_spectral.py``):
  ``compute_laplacian_spectrum_weigV(backend)``,
  ``compute_k_eigvV(k, backend, which)``,
  ``get_eigV(which)``, ``get_signed_laplacian_embedding(k)``

**Information theory** (``_infotheory.py``):
  ``compute_signed_laplacian_entropy(steps, t1, t2)``,
  ``get_entropy()``, ``get_specific_heat()``

**Clustering** (``_partitioning.py``):
  ``get_eigV_cluster_sizes(which)``, ``get_cluster_distribution(which)``,
  ``make_clustersYN(k, val)``, ``compute_pinf(which)``

**Dynamics/energy** (``_dynamics.py``):
  ``compute_rbim_energy_eigV(which)``, ``compute_sksph_energy_eigV(which)``,
  ``get_ferroAntiferro_regions(attr_str)``

**I/O** (``_exports.py`` / ``_loaders.py``):
  ``export_eigV_all()``, ``load_eigV_all()``, ``set_edgel_from_bin()``


Complete Graph Type Catalogue
-----------------------------

The table below lists every concrete graph type with its key parameters,
backend support, and facade import path.

.. list-table::
   :header-rows: 1
   :widths: 22 38 10 10 20

   * - Type
     - Key Parameters
     - NX
     - GT
     - Category
   * - ``SignedGraph``
     - ``G=``, ``pflip=``, ``seed=``
     - Yes
     - Yes
     - Base
   * - ``Lattice2D``
     - ``side1=``, ``geo=`` (sqr/tri/hex/oct_sqr), ``pbc=``
     - Yes
     - Yes
     - Lattice
   * - ``Lattice3D``
     - ``dim=``, ``geo=`` (sc/bcc/fcc), ``pbc=``
     - Yes
     - Yes
     - Lattice
   * - ``LatticeND``
     - ``shape=``, ``geo=`` (routes to 2D/3D), ``periodic=``
     - Yes
     - Yes
     - Lattice
   * - ``ErdosRenyi``
     - ``n=``, ``p=``, ``extract_giant_component=``
     - Yes
     - Yes
     - Random
   * - ``kRegularGraph``
     - ``n=``, ``k=``
     - Yes
     - FB
     - Random
   * - ``ConfigurationModel``
     - ``degree_sequence=``
     - Yes
     - FB
     - Random
   * - ``RandomGeometric``
     - ``n=``, ``radius=``, ``dim=``
     - Yes
     - FB
     - Random
   * - ``BarabasiAlbert``
     - ``n=``, ``m=``
     - Yes
     - FB
     - Scale-Free
   * - ``ExtendedBarabasiAlbert``
     - ``n=``, ``m=``, ``p=``, ``q=``
     - Yes
     - FB
     - Scale-Free
   * - ``DualBarabasiAlbert``
     - ``n=``, ``m1=``, ``m2=``, ``p=``
     - Yes
     - FB
     - Scale-Free
   * - ``HolmeKim``
     - ``n=``, ``m=``, ``p=``
     - Yes
     - FB
     - Scale-Free
   * - ``WattsStrogatz``
     - ``n=``, ``k=``, ``p=``
     - Yes
     - Yes
     - Small-World
   * - ``StochasticBlockModel``
     - ``sizes=``, ``p_matrix=``
     - Yes
     - Yes
     - Community
   * - ``LFRBenchmark``
     - ``n=``, ``tau1=``, ``tau2=``, ``mu=``
     - Yes
     - FB
     - Community
   * - ``FullyConnected``
     - ``N=`` (capital)
     - Yes
     - FB
     - Complete
   * - ``BipartiteGraph``
     - ``n1=``, ``n2=``, ``p=`` or ``top_degrees=``, ``bottom_degrees=``
     - Yes
     - Yes
     - Bipartite
   * - ``BipartiteFromDegreeSequence``
     - ``top_degrees=``, ``bottom_degrees=``
     - Yes
     - Yes
     - Bipartite
   * - ``HofieldNN``
     - ``with_patterns=``, ``digit=``, ``n_samples=``
     - Yes
     - FB
     - Neural
   * - ``SCSGeneralizedNN``
     - ``N=``, ``gamma=``, ``J0=``, ``J=``
     - Yes
     - FB
     - Neural
   * - ``MultiplicativeCascade``
     - ``p1=``, ``p2=``, ``fraction=``, ``iterations=``
     - Yes
     - Yes
     - Multispectral
   * - ``VicsekGraph``
     - ``N=``, ``k=``, ``pij=``
     - Yes
     - Yes
     - Multispectral
   * - ``HierarchicalModular``
     - ``levels=``, ``branching=``, ``leaf_nodes=``
     - Yes
     - FB
     - Multispectral
   * - ``GraphOfGraphs``
     - ``base_graph_type=``, ``fiber_graph_type=``
     - Yes
     - Yes
     - Graph-of-Graphs
   * - ``DiracCombGraph``
     - ``base_nodes=``, ``fiber_nodes=``
     - Yes
     - Yes
     - Graph-of-Graphs
   * - ``DiracBrushGraph``
     - ``base_x=``, ``base_y=``, ``fiber_nodes=``
     - Yes
     - Yes
     - Graph-of-Graphs
   * - ``SierpinskiGraph``
     - ``n=``
     - Yes
     - FB
     - Fractal
   * - ``DGMgraph``
     - ``n=``
     - Yes
     - FB
     - Fractal
   * - ``TemporalGraph``
     - ``n_nodes=``, ``n_snapshots=``, ``time_unit=``
     - Yes
     - FB
     - Temporal

**Legend:** Yes = native implementation, FB = fallback to NX with warning.

All types accept ``pflip=`` (float, 0.0–1.0), ``seed=`` (int), and
``engine=`` (str, ``'nx'``/``'gt'``).


LatticeND: Unified Lattice Dispatcher
--------------------------------------

``LatticeND`` is a smart entry point that routes to the appropriate
specialized implementation based on dimensionality and geometry:

.. code-block:: python

   from lrgsglib.graphs import LatticeND

   # 2D → dispatches to Lattice2D
   lat = LatticeND(shape=(64, 64), geo='tri', pflip=0.3)
   # type(lat) == Lattice2DNX

   # 3D → dispatches to Lattice3D
   lat = LatticeND(shape=(10, 10, 10), geo='bcc')
   # type(lat) == Lattice3DNX

   # 4D+ → generic cubic lattice
   lat = LatticeND(shape=(5, 5, 5, 5), periodic=True)
   # type(lat) == LatticeNDNX

   # Default geometry: sqr for 2D, sc for 3D
   lat = LatticeND(shape=(20, 20))  # same as Lattice2D(side1=20, geo='sqr')

The specialized classes (``Lattice2D``, ``Lattice3D``) remain available
for direct use when you need dimension-specific features like custom
geometries, positions, dilution, or cell defects.


TemporalGraph
-------------

``TemporalGraph`` represents time-varying signed networks. It inherits
from ``SignedGraphNX``, giving full access to spectral analysis, entropy,
and clustering on any materialized snapshot.

.. code-block:: python

   from lrgsglib.graphs import TemporalGraph

   tg = TemporalGraph(n_nodes=100, pflip=0.1, seed=42)

   # Add temporal edges
   tg.add_temporal_edge(0, 1, t_start=0, t_end=10, sign=1)
   tg.add_temporal_edge(1, 2, t_start=5, t_end=15, sign=-1)

   # Get plain snapshot (nx.Graph)
   G = tg.get_snapshot(time=7.0)

   # Get full SignedGraphNX snapshot (with spectral methods)
   sg = tg.get_signed_snapshot(time=7.0)
   sg.compute_laplacian_spectrum_weigV()
   entropy = sg.get_entropy()

   # Sign evolution over time
   tg.flip_random_edges(time=5.0, fraction=0.2)
   evolution = tg.get_frustration_evolution()

.. note::

   ``TemporalSignedGraphNX`` is deprecated. Use ``TemporalGraphNX``
   directly — it now inherits from ``SignedGraphNX`` and includes all
   sign-flipping methods.


Common Patterns
---------------

Creating a graph and computing spectrum:

.. code-block:: python

   from lrgsglib.graphs import ErdosRenyi

   er = ErdosRenyi(n=500, p=0.02, pflip=0.2, seed=42)
   er.flip_random_fract_edges()

   # Compute full spectrum
   er.compute_laplacian_spectrum_weigV()

   # Access results
   eigenvalues = er.eigv
   print(f"Spectral gap: {eigenvalues[1]:.4f}")

Accessing graph representations:

.. code-block:: python

   # Integer-labeled graph (always available)
   G = er.gr['G']

   # Coordinate-labeled graph (lattices only)
   from lrgsglib.graphs import Lattice2D
   lat = Lattice2D(side1=10, geo='sqr')
   H = lat.gr['H']  # nodes are (x, y) tuples

Switching engines:

.. code-block:: python

   from lrgsglib.graphs import BarabasiAlbert

   ba_nx = BarabasiAlbert(n=1000, m=3, engine='nx')
   ba_gt = BarabasiAlbert(n=1000, m=3, engine='gt')

   # Both satisfy SignedGraphProtocol
   print(ba_nx.N, ba_gt.N)  # same logical graph, different backends
