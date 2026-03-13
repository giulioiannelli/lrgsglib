Lattices
========

``lrgsglib`` provides 2D, 3D, and N-dimensional lattice generators that build
signed graphs with coordinate representations and plotting-ready positions.

- :py:class:`~lrgsglib.graphs.Lattice2D.Lattice2D` — 2D lattices with specialized geometries
- :py:class:`~lrgsglib.graphs.Lattice3D.Lattice3D` — 3D lattices with specialized geometries
- :py:class:`~lrgsglib.graphs.LatticeND.LatticeND` — Unified entry point (dispatches to 2D/3D or generic cubic)

2D Lattices
-----------

Create a square lattice and apply sign flips:

.. code-block:: python

   from lrgsglib.graphs import Lattice2D

   lat = Lattice2D(side1=20, geo="sqr", pflip=0.2, seed=1)
   lat.flip_random_fract_edges()

   print(lat.N, lat.Ne)

Common 2D geometries (short codes):

- ``"sqr"``: square lattice
- ``"tri"``: triangular lattice
- ``"hex"``: hexagonal lattice
- ``"sqr_sw"`` / ``"tri_sw"``: small-world variants
- ``"oct_sqr"``: rhomb-octagonal lattice

You can also pass the full geometry names (``"squared"``, ``"triangular"``,
``"hexagonal"``, ``"squared_sw"``, ``"triangular_sw"``, ``"octagonal_sqr"``).

3D Lattices
-----------

Create a simple cubic lattice:

.. code-block:: python

   from lrgsglib.graphs import Lattice3D

   lat3 = Lattice3D(dim=8, geo="sc", pflip=0.1, seed=3)
   lat3.flip_random_fract_edges()

   print(lat3.N, lat3.Ne)

Supported 3D geometries:

- ``"sc"`` or ``"simple_cubic"``
- ``"bcc"`` or ``"body_centered"``
- ``"fcc"`` or ``"face_centered"``

Graph Representations and Coordinates
-------------------------------------

Lattices expose two graph representations:

- ``H``: nodes labeled by their coordinate tuples
- ``G``: nodes relabeled to consecutive integers (default)

Most methods accept ``on_g=`` to select which representation to use.

Positions and Plotting
----------------------

Set ``with_positions=True`` to store node positions for plotting:

.. code-block:: python

   lat = Lattice2D(side1=10, geo="tri", with_positions=True)
   lat.flip_random_fract_edges()

For ``Lattice3D``, positions are projected to 2D using the ``theta`` and
``phi`` angles.

N-Dimensional Lattices
----------------------

For arbitrary dimensions or when you want a single entry point, use
``LatticeND``. It automatically routes to the specialized class when a
known geometry is requested:

.. code-block:: python

   from lrgsglib.graphs import LatticeND

   # 2D → dispatches to Lattice2D
   lat = LatticeND(shape=(64, 64), geo="tri")

   # 3D → dispatches to Lattice3D
   lat = LatticeND(shape=(10, 10, 10), geo="bcc")

   # 4D+ → generic cubic lattice
   lat = LatticeND(shape=(5, 5, 5, 5), periodic=True)

For dimensions beyond 3D, only cubic lattice geometry is available.

Next Steps
----------

- :doc:`graph_architecture` for the complete graph type catalogue
- :doc:`spectral` for Laplacian spectra and embeddings
- :doc:`plotting` for lattice visualization helpers
