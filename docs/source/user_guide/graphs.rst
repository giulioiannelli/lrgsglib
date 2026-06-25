Signed Graphs
=============

Signed graphs represent networks with positive and negative relationships.
In ``lrgsglib``, a signed graph is a NetworkX graph with a ``weight`` edge
attribute, and the :py:class:`lrgsglib.graphs.SignedGraph.SignedGraph`
class adds utility methods for sign flips, spectral analysis, and dynamics.

This section covers the core concepts and the most common creation paths.

Core Concepts
-------------

- **Edge signs**: positive weights represent cooperative links; negative
  weights represent antagonistic or frustrated links.
- **Signed edge sets**: ``SignedGraph`` tracks edge sets in ``eset``
  (all edges), ``fleset`` (negative edges), and ``lfeset`` (positive edges).
- **Disorder**: signs (and continuous couplings) are assigned via the
  :py:class:`~lrgsglib.graphs.Disorder` spec and **realized at construction**.
  ``flip_sel_edges`` still applies a hand-picked set of flips. See
  :doc:`disorder` for the full model.

Creating a SignedGraph
----------------------

From a NetworkX graph — ``pflip`` is applied during construction:

.. code-block:: python

   import networkx as nx
   from lrgsglib.graphs import SignedGraph

   G = nx.path_graph(6)
   sg = SignedGraph(G, pflip=0.25, seed=123)  # already signed

   print(sg.N, sg.Ne, sg.Ne_n)

From built-in generators:

.. code-block:: python

   from lrgsglib.graphs import ErdosRenyi

   er = ErdosRenyi(n=200, p=0.05, pflip=0.2, seed=7)  # already signed

   print(er.N, er.Ne, er.Ne_n)

.. note::

   Construction realizes the disorder, so a separate
   ``flip_random_fract_edges()`` call is no longer required (it remains an
   idempotent no-op for backward compatibility). Pass ``disorder=None`` to defer
   sign realization, or a :py:class:`~lrgsglib.graphs.Disorder` for structured /
   continuous couplings — see :doc:`disorder`.

Signed Graph Utilities
----------------------

Apply sign flips directly:

.. code-block:: python

   edges_to_flip = [(0, 1), (1, 2)]
   sg.flip_sel_edges(edges_to_flip)

Inspect positive/negative edge sets:

.. code-block:: python

   on_g = sg.on_g
   num_neg = len(sg.fleset[on_g])
   num_pos = len(sg.lfeset[on_g])
   print(num_neg, num_pos)

Graph Representations
---------------------

Some topologies expose multiple graph representations. Lattices, for example,
use ``H`` for coordinate labels and ``G`` for integer labels.

- The dictionary ``sg.gr`` stores all available representations.
- Most methods accept ``on_g=`` to pick the active representation.
- ``map_node`` and ``map_edge`` contain node/edge mappings between
  representations.

Reproducibility
---------------

Set a ``seed`` in constructors to make stochastic steps deterministic:

.. code-block:: python

   sg = SignedGraph(G, pflip=0.1, seed=42)  # disorder realized deterministically

Next Steps
----------

- :doc:`disorder` for the full disorder model (supports × coupling laws)
- :doc:`lattices` for 2D and 3D lattice construction
- :doc:`spectral` for Laplacian spectra and embeddings
- :doc:`dynamics` for statistical physics simulations
