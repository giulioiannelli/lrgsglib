Changelog
=========

All notable changes to lrgsglib will be documented in this file.

Version 0.1.0 (2025-XX-XX)
--------------------------

Initial release of lrgsglib.

**Features:**

- Signed graph generation (Erdős-Rényi, lattices, custom topologies)
- 2D lattice support (square, triangular, hexagonal geometries)
- 3D lattice support (cubic geometry)
- Spectral analysis tools (Laplacian spectrum, frustration measures)
- Contact process simulation (Python and C backends)
- Ising dynamics simulation (C backend)
- Voter model simulation (C backend)
- Visualization tools for signed graphs and lattices
- C/C++ extensions for high-performance simulations
- Comprehensive test suite
- Documentation with Sphinx

**Modules:**

- ``nx_patches``: NetworkX extensions for signed graphs
- ``utils``: Utility functions for spectral analysis and more
- ``statsys``: Statistical physics simulation systems
- ``plotlib``: Plotting and visualization
- ``config``: Configuration and program arguments

**Build System:**

- Makefile-based build with modular configuration
- Support for standalone and submodule usage
- Automatic environment configuration (.env generation)
- C/C++ compilation with pybind11

**Testing:**

- pytest-based test suite
- Quick and extended test scripts
- Performance benchmarks
- Doctest support in documentation

**Documentation:**

- Comprehensive Sphinx documentation
- User guide with tutorials
- API reference
- Developer guide
- Installation instructions

Development Roadmap
-------------------

Future versions will include:

**Planned Features:**

- Additional graph topologies
- More statistical physics models
- GPU acceleration support
- Parallel simulation capabilities
- Enhanced visualization tools
- Theory section with mathematical background
- Comprehensive Python code examples

**Under Consideration:**

- Python-based alternatives to all C simulators
- Dynamic network evolution
- Network reconstruction algorithms
- Integration with graph-tool for large-scale graphs
