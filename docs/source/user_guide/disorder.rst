Disorder: signs and couplings
=============================

Every signed graph carries *disorder* — the assignment of negative signs (or,
more generally, continuous couplings) to a subset of its edges. ``lrgsglib``
expresses disorder through a single engine-neutral spec, :py:class:`Disorder`,
that is **realized at construction** and travels on the graph object as
``sg.disorder``. The same spec behaves identically on the NetworkX and
graph-tool engines and on every graph family (lattices, Erdős–Rényi,
Barabási–Albert, Holme–Kim, Watts–Strogatz, …).

.. note::

   **Construction now applies the disorder.** ``SomeGraph(pflip=0.2)`` returns a
   graph that is *already signed* — you no longer need to call
   ``flip_random_fract_edges()`` afterwards. (That call still exists and is an
   idempotent no-op on an already-realized graph, so old scripts keep working.)
   To recover the old "select now, realize later" behaviour, pass
   ``disorder=None`` — see :ref:`disorder-deferring`.

The two axes
------------

A disorder realization factors into two orthogonal choices:

- **support** — *which* edges are affected.
- **coupling law** — *what value* those edges receive.

.. list-table::
   :header-rows: 1
   :widths: 50 18 32

   * - call
     - support
     - law
   * - ``Lattice2D(pflip=0.2)``
     - ``rand``
     - ``flip`` (negate a 20 % random fraction)
   * - ``BarabasiAlbert(pflip=0.05, disorder='randXERR')``
     - ``randXERR``
     - ``flip`` (negate star defects)
   * - ``Lattice2D(disorder=Disorder('all', law='gaussian', params=...))``
     - ``all``
     - ``gaussian`` (continuous SK couplings)

The ``Disorder`` spec
---------------------

.. code-block:: python

   from lrgsglib.graphs import Disorder

   Disorder(
       support='rand',     # nwDict key | 'all' | 'none'
       law='flip',         # 'flip' | 'gaussian' | 'uniform' | 'powerlaw' | callable
       pflip=0.0,          # fraction / intensity for random supports
       params={},          # keyword args forwarded to the coupling law
       seed=None,          # reserved; by default the graph RNG is used
   )

The ``disorder=`` constructor argument accepts three forms, coerced by
``as_disorder``:

- ``None`` — **defer**: select the random fraction but do not realize signs.
- a **string** — shorthand for ``Disorder(support=<str>, law='flip',
  pflip=<top-level pflip>)``. ``'none'`` / ``''`` are aliases for ``None``.
- a :py:class:`Disorder` instance — used as-is (its own ``pflip`` wins).

The default is ``disorder='rand'`` with the top-level ``pflip`` (default
``0.0``), so a plain ``SomeGraph()`` is unsigned and ``SomeGraph(pflip=0.2)``
flips a 20 % random fraction — fully backward compatible.

Supports
--------

The support vocabulary is the ``nwDict`` pattern set plus two specials:

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - support
     - edges affected
   * - ``rand``
     - a random ``pflip`` fraction (the classic disorder)
   * - ``all``
     - every edge
   * - ``none`` / ``None``
     - nothing (defer)
   * - ``randXERR``
     - star (XERR) defects seeded at random nodes
   * - ``randZERR``
     - elementary-cell (ZERR) defects seeded at random nodes
   * - ``single``
     - one central edge
   * - ``singleXERR``
     - the star around the central node
   * - ``singleZERR``
     - the elementary cell through the central node

The structured (``*XERR`` / ``*ZERR`` / ``single*``) supports route through the
graph's ``nwDict`` and are built on demand. **They now work on every graph
family, not only lattices.** Non-geometric graphs use a hub-based central edge
(highest-degree node) for the ``single*`` patterns; lattices use their
coordinate-aware centre. ``ZERR`` is the smallest cycle through a node — chosen
*canonically* (the lexicographically-smallest minimal cycle, so both engines pick
the same cell; oriented to a fixed +x/+y plaquette on the square lattice) — and it
degrades gracefully (no edges) on a tree/pendant.

.. code-block:: python

   from lrgsglib.graphs import BarabasiAlbert, Disorder

   # Star defects on a scale-free graph — previously lattice-only.
   ba = BarabasiAlbert(n=400, m=3, seed=7,
                       disorder=Disorder('randXERR', pflip=0.05))
   print('randXERR' in ba.nwDict, ba.Ne_n)

Combining disorders
-------------------

Several supports can be combined into one :py:class:`CompositeDisorder` via the
``disorder=`` argument. The headline is a **mixture** — split the ``pflip`` budget
between defect types (the classic "half-X / half-Z"):

.. code-block:: python

   from lrgsglib.graphs import Lattice2D, Disorder, CompositeDisorder

   # half the seed budget grows XERR stars, half grows ZERR cells (total ~pflip)
   g = Lattice2D(side1=64, geo='sqr', pflip=0.10, seed=7,
                 disorder={'randXERR': 0.5, 'randZERR': 0.5})

The combination *mode* follows from the shorthand or an explicit operator:

.. list-table::
   :header-rows: 1
   :widths: 34 22 44

   * - syntax
     - mode
     - meaning
   * - ``{'randXERR': .5, 'randZERR': .5}``
     - ``mixture``
     - partition the ``pflip`` seed budget by the weights
   * - ``[A, B]``
     - ``overlay``
     - apply each in turn (idempotent for ``flip``)
   * - ``A | B``
     - ``union``
     - flip the union of the two edge sets
   * - ``A & B``
     - ``intersection``
     - flip only edges in both
   * - ``A - B``
     - ``difference``
     - flip A's edges that are not B's
   * - ``A ^ B``
     - ``symmetric_difference``
     - flip edges in exactly one

.. code-block:: python

   A = Disorder('randXERR', pflip=0.05)
   B = Disorder('randZERR', pflip=0.05)
   g = Lattice2D(side1=64, geo='sqr', seed=7, disorder=(A ^ B))

A mixture can also blend laws — e.g. half the seeds flipped ``±1``, half given
continuous couplings:

.. code-block:: python

   disorder = [(Disorder('randXERR', law='flip'),                          0.5),
               (Disorder('randZERR', law='gaussian', params={'sigma': 1.5}), 0.5)]

The ``mixture`` planner draws its seeds from a fixed RNG, so it is identical on
NetworkX and graph-tool; the set-algebra modes reuse the ``nwDict`` random
patterns and inherit their per-engine seeding. ``mixture`` (v1) combines the
``randXERR`` / ``randZERR`` family; use ``overlay`` for other supports.

Coupling laws
-------------

``flip`` (the default) negates the sign, ``w → -w`` (discrete ±1 disorder).
For continuous couplings (Sherrington–Kirkpatrick / spin-glass style), pass a
distributional law and its ``params``:

==================  ===================================================  ============================
law                 draws                                                params (defaults)
==================  ===================================================  ============================
``gaussian``        i.i.d. ``N(mu, sigma)``                              ``mu=0.0``, ``sigma=1.0``
``uniform``         i.i.d. on ``[low, high]``                            ``low=-1.0``, ``high=1.0``
``powerlaw``        Pareto magnitude × random sign                       ``exponent=2.5``, ``signed=True``, ``xmin=1.0``
``callable``        ``fn(n, rng, **params) -> array`` of length ``n``    (your own)
==================  ===================================================  ============================

.. code-block:: python

   from lrgsglib.graphs import Lattice2D, Disorder

   # Continuous Gaussian couplings on every edge.
   sg = Lattice2D(side1=32, geo='sqr', seed=7,
                  disorder=Disorder('all', law='gaussian',
                                    params={'mu': 0.0, 'sigma': 1.0}))
   L = sg.get_signed_laplacian()   # symmetric, PSD with continuous couplings

The signed Laplacian uses absolute-value degrees, so continuous couplings need
no special handling downstream. On graph-tool a ``weight`` (double) edge
property is created only when a distributional law applies; the discrete ±1
``sign`` path is untouched.

Registering a custom law
~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the :py:func:`register_coupling` decorator to add a named law:

.. code-block:: python

   import numpy as np
   from lrgsglib.graphs import register_coupling, Disorder, ErdosRenyi

   @register_coupling('bimodal')
   def _bimodal(n, rng, j=1.0, p=0.5):
       """±j with probability p / (1-p)."""
       return np.where(rng.uniform(0, 1, n) < p, -j, j)

   er = ErdosRenyi(n=300, p=0.03, seed=7,
                   disorder=Disorder('all', law='bimodal', params={'j': 1.0, 'p': 0.3}))

Geometry primitives (writing a new support)
-------------------------------------------

The structured supports are assembled from a small set of engine-neutral,
coordinate-optional primitives in
:py:mod:`lrgsglib.graphs._shared._nw_geometry`. Reach for these to build a
disorder pattern outside the built-in vocabulary (a future ``register_support``
hook will let such patterns drop in without subclassing the container):

.. list-table::
   :header-rows: 1
   :widths: 46 54

   * - primitive
     - returns
   * - ``star_edges(sg, node, on_g)``
     - edges incident to ``node`` — the **XERR** star
   * - ``elementary_cell_edges(sg, node, on_g)``
     - the canonical smallest cycle through ``node`` — the **ZERR** cell
   * - ``minimal_cycles_through(sg, node, on_g)``
     - *every* minimal cycle through ``node``
   * - ``oriented_cell_edges(sg, node, on_g, posfn, box)``
     - the coordinate-oriented plaquette (a lattice's +x/+y face)
   * - ``hub_central_edge(sg, on_g)``
     - the highest-degree hub edge (generic centre)
   * - ``neighbors_at_distance(sg, node, R, on_g)`` / ``ball_edges(sg, center, R, on_g)``
     - the nodes / edges of the radius-``R`` ball around a node
   * - ``common_neighbors(sg, a, b, on_g)`` / ``induced_edges(sg, nodes, on_g)``
     - shared-neighbour set / induced-subgraph edges

Two convenience methods on the graph wrap this seam so the engines agree
node-for-node:

- ``sg.cell_edges(node, on_g)`` — the ZERR cell. The baseline is the deterministic
  lexicographically-smallest minimal cycle (engine-independent); the square lattice
  overrides it with a coordinate-oriented +x/+y plaquette, so distinct seed nodes
  own distinct cells in the bulk.
- ``sg.get_central_edge(on_g)`` — the centre for the ``single*`` patterns:
  coordinate-aware on lattices, the degree hub elsewhere.

.. code-block:: python

   from lrgsglib.graphs import Lattice2D
   from lrgsglib.graphs._shared._nw_geometry import star_edges

   sg = Lattice2D(side1=8, geo='sqr', seed=1)
   star_edges(sg, 27, 'G')     # the XERR pattern at node 27
   sg.cell_edges(27, 'G')      # the oriented ZERR cell at node 27

Registering a custom support
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:py:func:`register_support` adds a named support whose builder returns the edges
to act on. Registered supports resolve directly to edges (bypassing ``nwDict``)
and compose with any coupling law:

.. code-block:: python

   from lrgsglib.graphs import register_support, Disorder, BarabasiAlbert
   from lrgsglib.graphs._shared._nw_geometry import star_edges

   @register_support('hub_star')
   def _hub_star(sg, pflip, rng, on_g, k=1):
       """Flip the stars of the k highest-degree hubs."""
       deg = lambda n: len(list(sg.get_graph_neighbors(n, on_g)))
       hubs = sorted(sg.get_nodes_list(on_g), key=deg, reverse=True)[:k]
       return [e for h in hubs for e in star_edges(sg, h, on_g)]

   ba = BarabasiAlbert(n=400, m=3, seed=7,
                       disorder=Disorder('hub_star', support_params={'k': 3}))

Draw all randomness from ``rng`` and sort inputs before sampling so the two
engines agree for a fixed seed. ``support_params`` (distinct from the law's
``params``) is forwarded to the builder. A ``hubXERR`` support ships built-in as
a worked example.

.. _disorder-deferring:

Deferring realization
---------------------

Pass ``disorder=None`` to keep the old two-step workflow — the random fraction
is selected at construction but signs are realized only when you ask:

.. code-block:: python

   from lrgsglib.graphs import BarabasiAlbert

   ba = BarabasiAlbert(n=400, m=3, pflip=0.2, seed=7, disorder=None)
   print(ba.Ne_n)                 # 0 — nothing realized yet
   ba.flip_random_fract_edges()   # realize the deferred selection
   print(ba.Ne_n)                 # > 0

Reproducibility
---------------

Pass a ``seed`` to the constructor; the disorder realization (random support
selection and distributional draws) is deterministic for a fixed seed and
graph.
