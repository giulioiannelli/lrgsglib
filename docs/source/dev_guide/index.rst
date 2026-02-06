Developer Guide
===============

This guide is for developers who want to contribute to lrgsglib or understand its internal architecture.

.. toctree::
   :maxdepth: 2

   architecture
   build_system
   c_extensions
   testing
   contributing
   style_guide

Overview
--------

The developer guide covers:

**Architecture**
   System design, module organization, and design patterns used in lrgsglib

**Build System**
   Makefile structure, compilation process, and environment configuration

**C/C++ Extensions**
   Integration of C/C++ code, pybind11 bindings, and performance-critical components

**Testing**
   Test organization, running tests, writing new tests, and benchmarking

**Contributing**
   Pull request workflow, code review process, and release procedures

**Style Guide**
   Code style conventions, documentation requirements, and best practices

Getting Started with Development
---------------------------------

To set up a development environment:

.. code-block:: bash

   # Clone the repository
   git clone https://github.com/giulioiannelli/lrgsglib
   cd lrgsglib
   git submodule init && git submodule update

   # Create conda environment
   conda env create -f lrgsgenv.yml
   conda activate lrgsgenv

   # Build and install in editable mode with dev tools
   make all
   pip install -e ".[dev,docs]"

   # Run tests to verify
   pytest test/

Development Workflow
--------------------

1. Create a feature branch
2. Make changes following the style guide
3. Add tests for new functionality
4. Run the test suite and ensure all tests pass
5. Update documentation
6. Submit a pull request

Quick Reference
---------------

**Running tests:**

.. code-block:: bash

   pytest test/                    # Full test suite
   python test/quick_test.py       # Quick smoke test
   pytest test/test_contact_process.py  # Specific test file

**Code formatting:**

.. code-block:: bash

   black src/                      # Format Python code
   isort src/                      # Sort imports
   flake8 src/                     # Check style
   mypy src/                       # Type checking

**Building docs:**

.. code-block:: bash

   cd docs
   make html                       # Build HTML docs
   make doctest                    # Run doctest validation
   make livehtml                   # Live preview with auto-reload

**Rebuilding C extensions:**

.. code-block:: bash

   make clean
   make all

For detailed information, see the individual sections of this guide.
