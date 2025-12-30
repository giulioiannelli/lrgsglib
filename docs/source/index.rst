lrgsglib Documentation
======================

**lrgsglib** is a Python library implementing the theoretical tools of the **Laplacian Renormalization Group for Signed Graphs**. It provides utilities for building signed networks, running renormalization flows, and simulating statistical physics models such as the Ising model, contact process, and voter dynamics.

The library combines high-level Python interfaces with performance-critical C/C++ extensions for efficient simulations on large networks.

Project Scope
-------------

This project analyzes (signed) graph structure through the lens of the
(signed) Laplacian operator, with a particular focus on diffusion and
anomalous diffusion phenomena (including Anderson-localization-style
behavior). The central goal is to disentangle *structure* from *dynamics*:
we study how spectral information (eigenvalues/eigenvectors) of possibly
frustrated, disordered systems shapes and constrains dynamical processes.

To connect structure with dynamics, the library provides a family of
universality-class dynamics (e.g., ContactProcess, Ising, Voter models),
so users can simulate stochastic processes on the same underlying graphs
and compare how diffusion modes influence macroscopic behavior. For speed,
several dynamics have optimized C backends.

Project Layout
--------------

- ``src/lrgsglib`` contains the reusable Python package (graph objects,
  dynamics classes, utilities, plotting helpers, and bindings).
- ``src/`` contains standalone or compound programs used to compute or
  reproduce specific experiments and data products.
- ``src/lrgsglib/Ccore`` hosts performance-critical C/C++ extensions used
  by selected dynamics.

.. note::
   This documentation is for version |release|. For the latest development version,
   see the `GitHub repository <https://github.com/giulioiannelli/lrgsglib>`_.

Features
--------

- **Signed Graph Generation**: Create and manipulate signed networks (Erdős-Rényi, lattices, custom topologies)
- **Spectral Analysis**: Compute Laplacian spectra, frustration measures, and renormalization flows
- **Statistical Physics Simulations**: Efficient implementations of Ising dynamics, contact process, and voter models
- **Lattice Support**: 2D and 3D lattices with various boundary conditions
- **Visualization**: Specialized plotting functions for signed graphs and lattices
- **C/C++ Extensions**: High-performance simulators built with pybind11

Quick Links
-----------

- **Installation** - :doc:`installation` - Get started with installation instructions
- **Quick Start** - :doc:`quickstart` - Jump right in with a quick tutorial
- **User Guide** - :doc:`user_guide/index` - Comprehensive tutorials and examples
- **API Reference** - :doc:`api/index` - Detailed documentation of all modules

Table of Contents
-----------------

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   user_guide/index
   user_guide/graphs
   user_guide/lattices
   user_guide/spectral
   user_guide/dynamics
   user_guide/plotting
   user_guide/examples

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/index

.. toctree::
   :maxdepth: 2
   :caption: Developer Guide

   dev_guide/index
   dev_guide/architecture
   dev_guide/build_system
   dev_guide/c_extensions
   dev_guide/testing
   dev_guide/contributing
   dev_guide/style_guide

.. toctree::
   :maxdepth: 2
   :caption: Theory

   theory/index

.. toctree::
   :maxdepth: 1
   :caption: Additional Information

   changelog

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Citation
--------

If you use lrgsglib in your research, please cite:

.. code-block:: bibtex

   @software{lrgsglib,
     title = {lrgsglib: Laplacian Renormalization Group for Signed Graphs},
     author = {Iannelli, Giulio and Villegas, Pablo},
     year = {2025},
     url = {https://github.com/giulioiannelli/lrgsglib},
     version = {0.1.0}
   }

License
-------

This project is licensed under the MIT License - see the LICENSE file for details.

Authors
-------

- **Giulio Iannelli** - giulioiannelli.w@gmail.com
- **Pablo Villegas** - pablo.villegas@cref.it
