graphs - Unified Graph Interface
================================

The ``graphs`` module provides a unified API for creating signed graphs using
multiple backends (NetworkX, graph-tool, igraph).

.. currentmodule:: lrgsglib.graphs

.. note::

   This is the **recommended** module for new code. The ``nx_patches`` module
   is deprecated but still functional for backward compatibility.

Overview
--------

This module includes:

- **Engine Management**: Select between NetworkX, graph-tool, and igraph backends
- **Lattice Graphs**: 2D and 3D lattice generation with various geometries
- **Random Graphs**: Erdos-Renyi, Barabasi-Albert, Watts-Strogatz, SBM
- **Protocols**: Type-safe interfaces for signed graph implementations
- **Submodules**: Direct access to engine-specific implementations

Quick Start
-----------

.. code-block:: python

   from lrgsglib.graphs import Lattice2D, ErdosRenyi, set_default_engine

   # Create graphs with explicit engine
   lat_nx = Lattice2D(side1=100, geo='sqr', engine='nx')
   lat_gt = Lattice2D(side1=100, geo='sqr', engine='gt')

   # Or set global default
   set_default_engine('gt')
   lat = Lattice2D(side1=100, geo='sqr')  # Uses graph-tool

Engine Management
-----------------

Functions to configure and query the graph backend system.

.. autoclass:: GraphEngine
   :members:
   :undoc-members:

.. autofunction:: set_default_engine

.. autofunction:: get_default_engine

.. autofunction:: get_implementation

.. autofunction:: register_implementation

.. autofunction:: is_engine_available

.. autofunction:: available_engines

.. autofunction:: list_registered_types

.. autofunction:: list_engines_for_type

Protocols
---------

Type-safe interfaces that all graph implementations must satisfy.
These enable static type checking and runtime validation.

.. autoclass:: SignedGraphProtocol
   :members:

.. autoclass:: SpectralGraphProtocol
   :members:

.. autoclass:: LatticeGraphProtocol
   :members:

.. autofunction:: is_signed_graph

.. autofunction:: is_lattice_graph

Lattice Graphs
--------------

Lattice2D
~~~~~~~~~

2D lattice graphs with square, triangular, or hexagonal geometry.

.. autofunction:: Lattice2D

**Parameters:**

- ``side1`` (int): Size of first dimension
- ``side2`` (int, optional): Size of second dimension (defaults to side1)
- ``geo`` (str): Geometry type ('sqr', 'tri', 'hex')
- ``pflip`` (float): Fraction of edges to flip sign (0.0 to 1.0)
- ``periodic`` (bool): Whether lattice is periodic (default True)
- ``seed`` (int, optional): Random seed for reproducibility
- ``engine`` (str): Backend engine ('nx', 'gt', or None for default)

**Example:**

.. code-block:: python

   from lrgsglib.graphs import Lattice2D

   # Create 100x100 square lattice with 20% frustrated edges
   lat = Lattice2D(side1=100, geo='sqr', pflip=0.2)
   lat.flip_random_fract_edges()

   # Access the underlying graph
   G = lat.gr['G']  # Integer-labeled nodes
   H = lat.gr['H']  # Coordinate-labeled nodes

Lattice3D
~~~~~~~~~

3D lattice graphs with simple cubic, BCC, or FCC geometry.

.. autofunction:: Lattice3D

**Parameters:**

- ``dim`` (int or tuple): Dimensions of lattice
- ``geo`` (str): Geometry type ('cubic', 'bcc', 'fcc')
- ``pflip`` (float): Fraction of edges to flip sign
- ``periodic`` (bool): Whether lattice is periodic
- ``seed`` (int, optional): Random seed
- ``engine`` (str): Backend engine

Random Graph Generators
-----------------------

ErdosRenyi
~~~~~~~~~~

Erdos-Renyi random graphs with optional signed edges.

.. autofunction:: ErdosRenyi

**Parameters:**

- ``n`` (int): Number of nodes
- ``p`` (float): Edge probability
- ``pflip`` (float): Fraction of edges to flip sign
- ``extract_giant_component`` (bool): Extract largest connected component
- ``seed`` (int, optional): Random seed
- ``engine`` (str): Backend engine

**Example:**

.. code-block:: python

   from lrgsglib.graphs import ErdosRenyi

   # Create sparse random graph with 1000 nodes
   er = ErdosRenyi(n=1000, p=0.01, pflip=0.1)

BarabasiAlbert
~~~~~~~~~~~~~~

Barabasi-Albert scale-free networks using preferential attachment.

.. autofunction:: BarabasiAlbert

**Parameters:**

- ``n`` (int): Number of nodes
- ``m`` (int): Number of edges to attach from new node
- ``pflip`` (float): Fraction of edges to flip sign
- ``seed`` (int, optional): Random seed
- ``engine`` (str): Backend engine

WattsStrogatz
~~~~~~~~~~~~~

Watts-Strogatz small-world networks.

.. autofunction:: WattsStrogatz

**Parameters:**

- ``n`` (int): Number of nodes
- ``k`` (int): Each node connected to k nearest neighbors in ring
- ``p`` (float): Rewiring probability
- ``pflip`` (float): Fraction of edges to flip sign
- ``seed`` (int, optional): Random seed
- ``engine`` (str): Backend engine

StochasticBlockModel
~~~~~~~~~~~~~~~~~~~~

Stochastic Block Model with community structure.

.. autofunction:: StochasticBlockModel

**Parameters:**

- ``sizes`` (list): Sizes of each community
- ``p_matrix`` (array): Block probability matrix
- ``pflip`` (float): Fraction of edges to flip sign
- ``extract_giant_component`` (bool): Extract largest connected component
- ``seed`` (int, optional): Random seed
- ``engine`` (str): Backend engine

Submodules
----------

Direct Access to Engine-Specific Implementations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For advanced use cases requiring direct access to engine-specific classes:

**NetworkX implementations:**

.. code-block:: python

   from lrgsglib.graphs.nx import (
       SignedGraphNX,
       Lattice2DNX,
       Lattice3DNX,
       ErdosRenyiNX,
       BarabasiAlbertNX,
       WattsStrogatzNX,
       StochasticBlockModelNX,
   )

**graph-tool implementations:**

.. code-block:: python

   from lrgsglib.graphs.gt import (
       SignedGraphGT,
       Lattice2DGT,
       Lattice3DGT,
       ErdosRenyiGT,
       BarabasiAlbertGT,
       WattsStrogatzGT,
       StochasticBlockModelGT,
   )

.. note::

   Some graph-tool implementations may fall back to NetworkX if native
   implementations are not yet available. A warning is emitted in such cases.

NetworkX Module
~~~~~~~~~~~~~~~

.. automodule:: lrgsglib.graphs.nx
   :members:
   :undoc-members:
   :show-inheritance:

graph-tool Module
~~~~~~~~~~~~~~~~~

.. automodule:: lrgsglib.graphs.gt
   :members:
   :undoc-members:
   :show-inheritance:

Deprecation Notice
------------------

.. warning::

   The ``nx_patches`` module is **deprecated and will be removed in v2.0**.

   Update your imports:

   .. code-block:: python

      # OLD (will stop working)
      from lrgsglib.nx_patches import Lattice2D

      # NEW (use this)
      from lrgsglib.graphs import Lattice2D

All graph types previously in ``nx_patches`` are now available through:

- ``lrgsglib.graphs`` - Unified API with engine selection
- ``lrgsglib.graphs.nx`` - Direct NetworkX implementations
- ``lrgsglib.graphs.gt`` - Direct graph-tool implementations
