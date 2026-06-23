plotlib - Plotting and Visualization
=====================================

The ``plotlib`` module provides specialized plotting and visualization functions for signed graphs, lattices, and statistical physics simulations.

.. currentmodule:: lrgsglib.plotlib

Overview
--------

This module includes visualization tools for:

- **Lattice Plotting**: Visualization of 2D and 3D lattices with custom geometries
- **3D Plotting**: Advanced 3D visualizations for lattices and networks
- **Tilings**: Geometric tiling patterns and visualizations
- **Animations**: Animated visualizations of dynamics
- **Color Management**: Color maps, palettes, and scales for signed networks
- **Matplotlib Extensions**: Custom axes patches and formatters
- **LRG Visualizations**: Laplacian renormalization group specific plots

Lattice Visualization
---------------------

.. automodule:: lrgsglib.plotlib.lattices
   :members:
   :undoc-members:
   :show-inheritance:

3D Plotting
-----------

.. automodule:: lrgsglib.plotlib.plot3d
   :members:
   :undoc-members:
   :show-inheritance:

Tilings
-------

.. automodule:: lrgsglib.plotlib.tilings
   :members:
   :undoc-members:
   :show-inheritance:

Animation
---------

Engine- and graph-agnostic animation primitives (codec / output handling). This
is the bottom layer of the animation stack and knows nothing about graphs.

* **Structural, graph-type renderers** (e.g. the engine-agnostic 2D-lattice
  ``imshow`` movies ``animate_states`` / ``animate_largest_cluster``) live in
  :mod:`lrgsglib.graphs._shared.animation` and are bound as methods on the
  lattice classes -- call them as ``lat.animate_states(states)``.
* **Dynamics frame collection** (running a model to produce state snapshots)
  lives with the model, e.g. :mod:`lrgsglib.statsys.ContactProcess.frames`.

.. automodule:: lrgsglib.plotlib.animation
   :members:
   :undoc-members:
   :show-inheritance:

Color Management
----------------

Color Utilities
~~~~~~~~~~~~~~~

.. automodule:: lrgsglib.plotlib.color
   :members:
   :undoc-members:
   :show-inheritance:

Colormaps
~~~~~~~~~

.. automodule:: lrgsglib.plotlib.colormaps
   :members:
   :undoc-members:
   :show-inheritance:

Colorbars
~~~~~~~~~

.. automodule:: lrgsglib.plotlib.colorbars
   :members:
   :undoc-members:
   :show-inheritance:

Matplotlib Extensions
---------------------

Axes Patches
~~~~~~~~~~~~

.. automodule:: lrgsglib.plotlib.ax_patches
   :members:
   :undoc-members:
   :show-inheritance:

Formatters
~~~~~~~~~~

.. automodule:: lrgsglib.plotlib.formatter
   :members:
   :undoc-members:
   :show-inheritance:

LRG Plotting
------------

.. automodule:: lrgsglib.plotlib.lrg
   :members:
   :undoc-members:
   :show-inheritance:

Mathematical Plotting
---------------------

.. automodule:: lrgsglib.plotlib.mathplot
   :members:
   :undoc-members:
   :show-inheritance:

Chladni Patterns
----------------

.. automodule:: lrgsglib.plotlib.chladni
   :members:
   :undoc-members:
   :show-inheritance:

Constants
---------

.. automodule:: lrgsglib.plotlib.const_plotlib
   :members:
   :undoc-members:
   :show-inheritance:
