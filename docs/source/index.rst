lrgsglib documentation
=======================

**Version**: |release|

**Useful links**:
`Installation <installation.html>`_ |
`Source repository <https://github.com/giulioiannelli/lrgsglib>`_ |
`Issue tracker <https://github.com/giulioiannelli/lrgsglib/issues>`_

**lrgsglib** is a Python library implementing the theoretical tools of the
**Laplacian Renormalization Group for Signed Graphs**: building signed
networks, running renormalization flows, and simulating statistical-physics
models (Ising, contact process, voter dynamics) on top of them. It pairs
high-level Python interfaces with performance-critical C/C++ extensions for
efficient simulations on large networks.

.. grid:: 1 2 2 2
    :gutter: 3
    :class-container: sd-text-center sd-mt-4

    .. grid-item-card::
        :link: quickstart
        :link-type: doc
        :class-card: sd-rounded-3 sd-shadow-md

        :octicon:`rocket;3em;sd-text-primary`

        **Getting started**
        ^^^
        New to lrgsglib? Install it, build your first signed graph and run a
        simulation in a handful of lines.

    .. grid-item-card::
        :link: user_guide/index
        :link-type: doc
        :class-card: sd-rounded-3 sd-shadow-md

        :octicon:`book;3em;sd-text-primary`

        **User guide**
        ^^^
        In-depth tutorials on graphs, lattices, spectral analysis, dynamics
        and plotting, with worked examples.

    .. grid-item-card::
        :link: api/index
        :link-type: doc
        :class-card: sd-rounded-3 sd-shadow-md

        :octicon:`code-square;3em;sd-text-primary`

        **API reference**
        ^^^
        Detailed description of every module, class and function: graphs,
        statsys dynamics, spectral utilities and bindings.

    .. grid-item-card::
        :link: dev_guide/index
        :link-type: doc
        :class-card: sd-rounded-3 sd-shadow-md

        :octicon:`tools;3em;sd-text-primary`

        **Developer guide**
        ^^^
        Architecture, the build system, C/C++ extensions, testing and the
        contribution conventions.

Project scope
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

Features
--------

- **Signed graph generation**: create and manipulate signed networks (Erdős-Rényi, lattices, custom topologies)
- **Spectral analysis**: compute Laplacian spectra, frustration measures, and renormalization flows
- **Statistical-physics simulations**: efficient implementations of Ising dynamics, contact process, and voter models
- **Lattice support**: 2D and 3D lattices with various boundary conditions
- **Visualization**: specialized plotting functions for signed graphs and lattices
- **C/C++ extensions**: high-performance simulators built with pybind11

.. Hidden top-level toctrees: they populate the navbar dropdowns and the
   (expandable) left sidebar without printing a long list on this page. Each
   section index owns the deeper pages, so nav expands one level at a time.

.. toctree::
   :hidden:
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :hidden:
   :caption: User Guide

   user_guide/index

.. toctree::
   :hidden:
   :caption: API Reference

   api/index

.. toctree::
   :hidden:
   :caption: Developer Guide

   dev_guide/index

.. toctree::
   :hidden:
   :caption: Theory

   theory/index

.. toctree::
   :hidden:
   :caption: Additional Information

   changelog

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
