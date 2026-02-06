Contributing
============

This document describes how to contribute to lrgsglib, including the
development workflow, pull request process, and release procedures.

Getting Started
---------------

Development Setup
~~~~~~~~~~~~~~~~~

1. **Fork and clone** the repository:

   .. code-block:: bash

      git clone https://github.com/YOUR_USERNAME/lrgsglib
      cd lrgsglib
      git remote add upstream https://github.com/giulioiannelli/lrgsglib

2. **Initialize submodules:**

   .. code-block:: bash

      git submodule init
      git submodule update

3. **Create conda environment:**

   .. code-block:: bash

      conda env create -f lrgsgenv.yml
      conda activate lrgsgenv

4. **Build and install:**

   .. code-block:: bash

      make all
      pip install -e ".[dev,docs]"

5. **Verify setup:**

   .. code-block:: bash

      pytest test/

Development Workflow
--------------------

Branch Strategy
~~~~~~~~~~~~~~~

- ``main`` - Stable release branch
- ``dev`` - Development branch (default for PRs)
- ``feature/description`` - Feature branches
- ``fix/description`` - Bug fix branches
- ``docs/description`` - Documentation branches

Creating a Feature Branch
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   # Update dev branch
   git checkout dev
   git pull upstream dev

   # Create feature branch
   git checkout -b feature/my-new-feature

   # Make changes...
   # Commit changes...

   # Push to your fork
   git push origin feature/my-new-feature

Making Changes
--------------

Code Changes
~~~~~~~~~~~~

1. **Follow the style guide** - See :doc:`style_guide`
2. **Add type hints** - All new code must have type annotations
3. **Write docstrings** - NumPy format with examples
4. **Add tests** - New functionality requires tests

Example commit message:

.. code-block:: text

   Add spectral entropy computation for signed graphs

   - Implement compute_spectral_entropy() in utils/lrg/infocomm.py
   - Add unit tests in test/test_spectral_entropy.py
   - Update user guide with examples

   Fixes #42

Documentation Changes
~~~~~~~~~~~~~~~~~~~~~

Documentation lives in two places:

1. **Library documentation** - ``docs/source/*.rst``
2. **Docstrings** - In Python source files

When adding features:

1. Update relevant ``.rst`` files
2. Add docstrings with examples
3. Run ``make -C docs html`` to verify

C/C++ Changes
~~~~~~~~~~~~~

1. Follow existing patterns in ``Ccore/``
2. Add Makefile rules for new binaries
3. Create Python wrappers in ``statsys/``
4. Benchmark before and after

Pull Request Process
--------------------

Before Submitting
~~~~~~~~~~~~~~~~~

1. **Run tests:**

   .. code-block:: bash

      pytest test/

2. **Check formatting:**

   .. code-block:: bash

      black --check src/
      isort --check src/
      flake8 src/
      mypy src/

3. **Build documentation:**

   .. code-block:: bash

      make -C docs html

4. **Update CHANGELOG** if applicable

Submitting a PR
~~~~~~~~~~~~~~~

1. Push your branch to your fork
2. Open a pull request against ``dev`` branch
3. Fill in the PR template:

   - **Summary** - What does this PR do?
   - **Changes** - List of changes
   - **Testing** - How was it tested?
   - **Related issues** - Link to issues

PR Review
~~~~~~~~~

- Address reviewer feedback
- Keep commits clean (squash if needed)
- Ensure CI passes

After Merge
~~~~~~~~~~~

.. code-block:: bash

   # Update your local dev
   git checkout dev
   git pull upstream dev

   # Delete feature branch
   git branch -d feature/my-new-feature

Code Review Guidelines
----------------------

For Reviewers
~~~~~~~~~~~~~

- Check code style and conventions
- Verify tests cover new functionality
- Ensure documentation is updated
- Consider performance implications
- Look for edge cases

For Authors
~~~~~~~~~~~

- Respond to all comments
- Explain design decisions
- Be open to suggestions
- Keep PRs focused and small

Issue Guidelines
----------------

Reporting Bugs
~~~~~~~~~~~~~~

Include:

1. **Description** - What happened?
2. **Expected behavior** - What should happen?
3. **Steps to reproduce** - Minimal example
4. **Environment** - OS, Python version, etc.
5. **Error messages** - Full traceback

Feature Requests
~~~~~~~~~~~~~~~~

Include:

1. **Use case** - Why is this needed?
2. **Proposed solution** - How might it work?
3. **Alternatives** - Other approaches considered
4. **Priority** - How important is this?

Release Process
---------------

Version Numbering
~~~~~~~~~~~~~~~~~

lrgsglib uses semantic versioning: ``MAJOR.MINOR.PATCH``

- **MAJOR** - Breaking changes
- **MINOR** - New features (backward compatible)
- **PATCH** - Bug fixes

Release Checklist
~~~~~~~~~~~~~~~~~

1. Update version in ``pyproject.toml``
2. Update ``CHANGELOG.rst``
3. Create release branch from ``dev``
4. Run full test suite
5. Build and verify documentation
6. Merge to ``main``
7. Tag release
8. Push to PyPI (if applicable)

Project Structure for Contributors
----------------------------------

Key directories to understand:

.. code-block:: text

   lrgsglib/
   ├── src/lrgsglib/          ← Main library code
   │   ├── nx_patches/        ← Graph classes (start here for new graphs)
   │   ├── statsys/           ← Dynamics (start here for new simulations)
   │   ├── utils/             ← Utility functions
   │   └── Ccore/             ← C code (performance-critical)
   │
   ├── src/                   ← Standalone programs
   │   ├── kernels/           ← Reusable computation blocks
   │   └── parsers/           ← Argument parsing
   │
   ├── test/                  ← Tests and benchmarks
   ├── docs/                  ← Sphinx documentation
   └── build/                 ← Build configuration

Common Contribution Scenarios
-----------------------------

Adding a New Graph Type
~~~~~~~~~~~~~~~~~~~~~~~

1. Create ``nx_patches/NewGraph/NewGraph.py``
2. Inherit from ``SignedGraph``
3. Add ``__init__.py`` with exports
4. Create ``kernels/NewGraph.py``
5. Add tests in ``test/test_new_graph.py``
6. Document in ``docs/source/api/``

Adding a New Simulation
~~~~~~~~~~~~~~~~~~~~~~~

1. Create ``statsys/NewDynamics.py``
2. Inherit from ``BinDynSys``
3. Implement Python version
4. Add C version in ``Ccore/statsys/``
5. Add Makefile rules
6. Add tests and benchmarks

Improving Performance
~~~~~~~~~~~~~~~~~~~~~

1. Profile to find bottleneck
2. Create benchmark in ``test/bench_*.py``
3. Implement optimization
4. Verify correctness with tests
5. Compare benchmark results

Getting Help
------------

- **GitHub Issues** - For bugs and features
- **Documentation** - This guide and API docs
- **Code Comments** - Check inline documentation

Questions about contributing? Open an issue with the "question" label.

See Also
--------

- :doc:`style_guide` - Code conventions
- :doc:`testing` - Writing tests
- :doc:`architecture` - System design
