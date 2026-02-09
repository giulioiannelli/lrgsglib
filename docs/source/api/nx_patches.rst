nx_patches (Deprecated)
=======================

.. deprecated:: 1.0

   The ``nx_patches`` module is deprecated. Use :mod:`lrgsglib.graphs.nx` instead.

All graph types previously available through ``nx_patches`` have been moved to
``lrgsglib.graphs.nx`` with the NX suffix naming convention.

Migration Guide
---------------

.. code-block:: python

   # Old (deprecated) - still works but emits DeprecationWarning
   from lrgsglib.nx_patches import Lattice2D, ErdosRenyi, SignedGraph

   # New (recommended) - use NX suffix classes
   from lrgsglib.graphs.nx import Lattice2DNX, ErdosRenyiNX, SignedGraphNX

   # Backward-compatible short names also available:
   from lrgsglib.graphs.nx import Lattice2D, ErdosRenyi, SignedGraph

   # Or use the unified multi-engine factory:
   from lrgsglib.graphs import Lattice2D
   lat = Lattice2D(10, engine='nx')  # NetworkX backend
   lat = Lattice2D(10, engine='gt')  # graph-tool backend

See :doc:`graphs` for the full API reference of the new module.
