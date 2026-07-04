graphs - Unified Graph Interface
================================

The ``graphs`` module provides a unified API for creating signed graphs using
multiple backends (NetworkX, graph-tool).

.. currentmodule:: lrgsglib.graphs

.. note::

   All facades are **classes** whose ``__new__`` dispatches to the correct
   engine implementation. They accept ``engine='nx'`` or ``engine='gt'``.

Overview
--------

This module includes:

- **Engine Management**: Select between NetworkX and graph-tool backends
- **Facade Classes**: Engine-agnostic entry points for all graph types
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

.. autoclass:: SignedGraphProtocol
   :members:

.. autoclass:: SpectralGraphProtocol
   :members:

.. autoclass:: LatticeGraphProtocol
   :members:

.. autofunction:: is_signed_graph

.. autofunction:: is_lattice_graph


Base
----

.. autoclass:: SignedGraph

Disorder
--------

Edge signs and continuous couplings are assigned through the engine-neutral
``Disorder`` spec (support × coupling-law), realized at construction. See the
:doc:`../user_guide/disorder` guide for usage.

.. autoclass:: Disorder
   :members:

.. autofunction:: register_coupling

Lattice Graphs
--------------

.. autoclass:: Lattice2D

.. autoclass:: Lattice3D

.. autoclass:: LatticeND

Random Graph Generators
-----------------------

.. autoclass:: ErdosRenyi

.. autoclass:: kRegularGraph

.. autoclass:: ConfigurationModel

.. autoclass:: RandomGeometric

Scale-Free Graphs
-----------------

.. autoclass:: BarabasiAlbert

.. autoclass:: ExtendedBarabasiAlbert

.. autoclass:: DualBarabasiAlbert

.. autoclass:: HolmeKim

Small-World Graphs
------------------

.. autoclass:: WattsStrogatz

Community Structure
-------------------

.. autoclass:: StochasticBlockModel

.. autoclass:: LFRBenchmark

Complete Graphs
---------------

.. autoclass:: FullyConnected

Bipartite Graphs
----------------

.. autoclass:: BipartiteGraph

.. autoclass:: BipartiteFromDegreeSequence

Neural Network Graphs
---------------------

.. autoclass:: HofieldNN

.. autoclass:: SCSGeneralizedNN

Multispectral Graphs
--------------------

.. autoclass:: MultiplicativeCascade

.. autoclass:: VicsekGraph

.. autoclass:: HierarchicalModular

Fractal Graphs
--------------

.. autoclass:: SierpinskiGraph

.. autoclass:: DGMgraph

Graph-of-Graphs
---------------

.. autoclass:: GraphOfGraphs

.. autoclass:: DiracCombGraph

.. autoclass:: DiracBrushGraph

Temporal Graphs
---------------

.. autoclass:: TemporalGraph


Submodules
----------

Direct Access to Engine-Specific Implementations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**NetworkX implementations:**

.. code-block:: python

   from lrgsglib.graphs.nx import (
       SignedGraphNX, LatticeNDNX,
       Lattice2DNX, Lattice3DNX,
       ErdosRenyiNX, BarabasiAlbertNX,
       WattsStrogatzNX, StochasticBlockModelNX,
       GraphOfGraphsNX, TemporalGraphNX,
   )

**graph-tool implementations:**

.. code-block:: python

   from lrgsglib.graphs.gt import (
       SignedGraphGT, LatticeNDGT,
       Lattice2DGT, Lattice3DGT,
       ErdosRenyiGT, BarabasiAlbertGT,
       WattsStrogatzGT, StochasticBlockModelGT,
       GraphOfGraphsGT,
   )

.. note::

   Some graph-tool implementations fall back to NetworkX if native
   implementations are not yet available. A warning is emitted in such cases.

NetworkX Module
~~~~~~~~~~~~~~~

.. automodule:: lrgsglib.graphs.nx
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

graph-tool Module
~~~~~~~~~~~~~~~~~

.. automodule:: lrgsglib.graphs.gt
   :members:
   :undoc-members:
   :show-inheritance:
   :no-index:

Module Structure
----------------

- ``lrgsglib.graphs`` - Unified API with engine selection
- ``lrgsglib.graphs.nx`` - Direct NetworkX implementations
- ``lrgsglib.graphs.gt`` - Direct graph-tool implementations
