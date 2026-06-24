API Reference
=============

This section contains the detailed API documentation for all modules, classes, and functions in lrgsglib.

.. toctree::
   :maxdepth: 2

   graphs
   utils
   statsys
   plotlib
   bindings
   config

Module Overview
---------------

**graphs**
   Unified graph interface with multi-engine support (NetworkX, graph-tool).
   Provides lattice graphs, random graphs, and signed graph operations.

**utils**
   Utility functions for spectral analysis, basic operations, and reconstruction algorithms

**statsys**
   Statistical physics simulation systems (Ising, contact process, voter models)

**plotlib**
   Plotting and visualization functions specialized for signed graphs and
   lattices, including the engine- and graph-agnostic animation primitives
   (``plotlib.animation``). Structural, graph-type animation renderers live
   under ``graphs._shared.animation`` (exposed via the ``lat.animate`` /
   ``lat.plot`` accessors on the graph classes) and dynamics frame collection
   under ``statsys`` (e.g. ``statsys.ContactProcess.frames``).

**bindings**
   Native language bindings for performance-critical operations

**config**
   Configuration management, error handling, and program arguments

Quick Links
-----------

Common modules:

- :py:mod:`lrgsglib.graphs` - Unified graph interface (Lattice2D, ErdosRenyi, etc.)
- :py:mod:`lrgsglib.graphs.nx` - NetworkX implementations
- :py:mod:`lrgsglib.graphs.gt` - graph-tool implementations
- :py:mod:`lrgsglib.utils.lrg.spectral` - Spectral analysis functions
- :py:mod:`lrgsglib.statsys` - Statistical physics simulations

Index
-----

* :ref:`genindex`
* :ref:`modindex`
