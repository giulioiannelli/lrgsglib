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
coordinate-aware centre. ``ZERR`` is the smallest cycle through a node, so it
degrades gracefully (no edges) on a tree/pendant.

.. code-block:: python

   from lrgsglib.graphs import BarabasiAlbert, Disorder

   # Star defects on a scale-free graph — previously lattice-only.
   ba = BarabasiAlbert(n=400, m=3, seed=7,
                       disorder=Disorder('randXERR', pflip=0.05))
   print('randXERR' in ba.nwDict, ba.Ne_n)

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
