# Documentation Build Plan for lrgsglib

**Project**: lrgsglib - Laplacian Renormalization Group for Signed Graphs
**Plan Created**: 2025-12-08
**Status**: Pending Approval

---

## Executive Summary

Build a comprehensive, reproducible, and testable documentation system for the lrgsglib library using Sphinx. The documentation will serve both research users and developers, with automated testing, versioning, and deployment capabilities.

---

## Current State Assessment

### Existing Documentation
- **README.md**: Basic installation and usage instructions
- **CLAUDE.md**: AI assistant guidance (project structure, conventions)
- **AGENTS.md**: Development conventions and code style
- **Scattered READMEs**: In `Ccore/` subdirectories (RBIsingM, voterM, signedRw, etc.)
- **Docstrings**: Mixed quality - some functions have NumPy-style docstrings, many lack documentation
- **No formal documentation system**: No Sphinx, no API reference, no user guides

### Current Tooling
- **Testing**: pytest-based (working, comprehensive tests in `test/`)
- **Code quality**: black, isort, flake8, mypy configured in pyproject.toml
- **Dependencies**: Listed in pyproject.toml with optional dev/jupyter extras
- **Build system**: Custom Makefile with modular configuration

### Documentation Needs
1. **User Documentation**: Installation, tutorials, examples, API reference
2. **Developer Documentation**: Architecture, C/C++ integration, contributing guide
3. **Reproducibility**: Doctests, example validation, automated builds
4. **Versioning**: Track API changes across versions
5. **Accessibility**: Searchable, browsable, hosted (GitHub Pages or ReadTheDocs)

---

## Documentation Architecture

### Technology Stack

**Core Tool: Sphinx**
- Industry standard for Python documentation
- Supports NumPy/Google docstring styles
- Extensible with plugins
- Integrates with testing frameworks

**Key Extensions:**
- `sphinx.ext.autodoc` - Auto-generate API docs from docstrings
- `sphinx.ext.autosummary` - Generate summary tables
- `sphinx.ext.doctest` - Test code examples in documentation
- `sphinx.ext.napoleon` - Parse NumPy/Google style docstrings
- `sphinx.ext.intersphinx` - Link to external docs (NumPy, NetworkX, etc.)
- `sphinx.ext.viewcode` - Add source code links
- `sphinx.ext.mathjax` - Render mathematical equations
- `sphinx_rtd_theme` - ReadTheDocs theme (clean, searchable)
- `myst-parser` - Support Markdown alongside reStructuredText
- `nbsphinx` - Include Jupyter notebooks as documentation

**Testing Integration:**
- Doctest: validate code examples in docstrings
- pytest-doctestplus: enhanced doctest support
- Continuous validation of examples

**Version Control:**
- sphinx-multiversion: support multiple doc versions
- Git-based versioning aligned with releases

---

## Documentation Structure

```
lrgsglib/
├── docs/                          # NEW: Documentation root
│   ├── source/                    # Sphinx source files
│   │   ├── conf.py               # Sphinx configuration
│   │   ├── index.rst             # Documentation home
│   │   ├── installation.rst      # Installation guide
│   │   ├── quickstart.rst        # Quick start tutorial
│   │   ├── user_guide/           # User documentation
│   │   │   ├── index.rst
│   │   │   ├── graphs.rst        # Signed graph basics
│   │   │   ├── lattices.rst      # Lattice generation
│   │   │   ├── spectral.rst      # Spectral analysis
│   │   │   ├── dynamics.rst      # Statistical physics simulations
│   │   │   ├── plotting.rst      # Visualization
│   │   │   └── examples.rst      # Complete examples
│   │   ├── api/                  # API reference (auto-generated)
│   │   │   ├── index.rst
│   │   │   ├── nx_patches.rst    # NetworkX extensions
│   │   │   ├── utils.rst         # Utilities
│   │   │   ├── statsys.rst       # Statistical systems
│   │   │   ├── plotlib.rst       # Plotting library
│   │   │   └── config.rst        # Configuration
│   │   ├── dev_guide/            # Developer documentation
│   │   │   ├── index.rst
│   │   │   ├── architecture.rst  # System architecture
│   │   │   ├── build_system.rst  # Build process
│   │   │   ├── c_extensions.rst  # C/C++ integration
│   │   │   ├── testing.rst       # Testing guide
│   │   │   ├── contributing.rst  # Contribution guidelines
│   │   │   └── style_guide.rst   # Code style (from AGENTS.md)
│   │   ├── notebooks/            # Tutorial notebooks
│   │   │   └── *.ipynb           # Jupyter notebook examples
│   │   ├── theory/               # Theoretical background
│   │   │   ├── index.rst
│   │   │   ├── lrg_theory.rst    # LRG mathematical framework
│   │   │   ├── signed_graphs.rst # Signed graph theory
│   │   │   └── stat_phys.rst     # Statistical physics models
│   │   ├── changelog.rst         # Version history
│   │   └── _static/              # Static assets (images, CSS)
│   ├── build/                    # Generated HTML/PDF (gitignored)
│   ├── Makefile                  # Sphinx build commands
│   └── make.bat                  # Windows build script
├── pyproject.toml                # UPDATED: Add docs dependencies
└── README.md                     # UPDATED: Link to full docs
```

---

## Implementation Plan

### Phase 1: Setup and Infrastructure (Reproducible Foundation)

**Goal**: Establish a working Sphinx documentation system with testing

**Tasks**:

1. **Install and configure Sphinx**
   - Add documentation dependencies to pyproject.toml
   - Create docs/ directory structure
   - Initialize Sphinx with `sphinx-quickstart`
   - Configure conf.py with extensions and theme

2. **Set up doctest integration**
   - Enable `sphinx.ext.doctest` extension
   - Configure pytest to run doctests
   - Create doctest conventions document
   - Add `make doctest` target

3. **Configure build automation**
   - Create docs/Makefile with targets: html, pdf, doctest, clean, livehtml
   - Add GitHub Actions workflow for docs CI/CD
   - Set up ReadTheDocs or GitHub Pages integration
   - Add docs build to main test suite

4. **Create base documentation files**
   - Write initial index.rst with project overview
   - Write installation.rst (migrate from README)
   - Create quickstart.rst with minimal example
   - Set up API structure (empty for now)

**Deliverables**:
- Working Sphinx build system
- Automated doctest validation
- CI/CD pipeline for docs
- Initial skeleton documentation

**Testing Criteria**:
- `make html` builds without errors
- `make doctest` runs successfully (even if no tests yet)
- Documentation displays correctly locally
- CI/CD deploys to hosting platform

---

### Phase 2: API Documentation (Auto-generated Reference)

**Goal**: Comprehensive API reference with validated examples

**Tasks**:

1. **Set up autodoc for Python modules**
   - Configure autodoc options in conf.py
   - Create API reference structure (api/ directory)
   - Generate autosummary tables for each module
   - Add module-level documentation

2. **Document core modules with examples**
   - `nx_patches/`: Lattice2D, Lattice3D, SignedGraph, ErdosRenyi
   - `utils/lrg/`: spectral, clustering
   - `utils/basic/`: common utilities
   - `statsys/`: ContactProcess, IsingDynamics, VoterModel
   - `plotlib/`: visualization functions

3. **Add and validate docstring examples**
   - Review existing docstrings (many already have NumPy style)
   - Add missing docstrings to public functions
   - Add executable examples to all major functions
   - Run doctest to validate examples

4. **Document C/C++ integration**
   - Create dev_guide/c_extensions.rst
   - Document pybind11 bindings
   - Explain C program interfaces
   - Document SFMT RNG usage

**Deliverables**:
- Complete API reference for all public modules
- Validated docstring examples
- C/C++ integration documentation

**Testing Criteria**:
- All public functions have docstrings
- All docstring examples pass doctest
- API docs render correctly with cross-references
- Links to source code work

---

### Phase 3: User Guide (Tutorials and Examples)

**Goal**: Comprehensive user-facing documentation

**Tasks**:

1. **Write conceptual guides**
   - Signed graphs: creation, properties, visualization
   - Lattices: 2D/3D generation, boundary conditions
   - Spectral analysis: Laplacian spectrum, renormalization
   - Statistical physics: Ising, contact process, voter model
   - Plotting: graph visualization, heatmaps, animations

2. **Create end-to-end examples**
   - Example 1: Creating and analyzing a signed Erdős-Rényi graph
   - Example 2: Running contact process on a 2D lattice
   - Example 3: Computing Laplacian spectrum and frustration
   - Example 4: Visualizing Ising dynamics with animations
   - Example 5: Submodule usage in lrgsglib-ipynb context

3. **Integrate Jupyter notebooks**
   - Select representative notebooks from ipynb/devTool/
   - Clean and annotate notebooks for documentation
   - Use nbsphinx to include in docs
   - Add narrative explanations

4. **Write theory section**
   - Laplacian renormalization group overview
   - Signed graph theory basics
   - Statistical physics model descriptions
   - References to academic papers

**Deliverables**:
- Complete user guide with tutorials
- 5+ end-to-end examples with validated code
- Integrated Jupyter notebooks
- Theory section with mathematical background

**Testing Criteria**:
- All tutorial code executes successfully
- Examples are reproducible with specified seeds
- Notebooks render correctly in docs
- Theory equations display properly

---

### Phase 4: Developer Documentation

**Goal**: Enable contributors to understand and extend the library

**Tasks**:

1. **Architecture documentation**
   - System overview diagram
   - Module dependency graph
   - Data flow in simulations
   - File organization rationale

2. **Build system documentation**
   - Makefile targets and configuration
   - Environment variable setup (.env generation)
   - Submodule vs standalone builds
   - C/C++ compilation process

3. **Testing guide**
   - Test organization (unit, integration, performance)
   - Running tests (pytest, quick_test, extended_test)
   - Writing new tests
   - Benchmark interpretation

4. **Contributing guide**
   - Code style enforcement (black, isort, mypy, flake8)
   - Pull request workflow
   - Documentation requirements
   - Release process

**Deliverables**:
- Complete developer guide
- Build system documentation
- Testing and contribution guidelines
- Architecture diagrams

**Testing Criteria**:
- New contributors can build from source following docs
- Code style guide is clear and enforced
- Test documentation covers all test types

---

### Phase 5: Polish and Deployment

**Goal**: Production-ready documentation system

**Tasks**:

1. **Documentation testing and validation**
   - Run complete doctest suite
   - Check all cross-references
   - Validate external links (intersphinx)
   - Test on fresh Python environment

2. **Accessibility and usability**
   - Add search functionality tuning
   - Create comprehensive index
   - Add glossary of terms
   - Improve navigation (TOC depth, sidebar)

3. **Versioning setup**
   - Configure sphinx-multiversion
   - Tag version 0.1.0 for initial release
   - Set up version switcher in theme
   - Document versioning policy

4. **Deployment automation**
   - Finalize GitHub Pages / ReadTheDocs setup
   - Add deployment to release workflow
   - Create documentation update guide
   - Set up automated link checking

5. **Migration and cleanup**
   - Update main README.md to point to docs
   - Archive or integrate scattered READMEs
   - Add "Documentation" badge to README
   - Announce documentation in project

**Deliverables**:
- Fully tested documentation build
- Versioned documentation with switcher
- Automated deployment pipeline
- Updated project README

**Testing Criteria**:
- Documentation builds with zero warnings
- All doctests pass
- Deployment succeeds automatically
- Version switcher works correctly

---

## Documentation Standards and Conventions

### Docstring Style

**Standard**: NumPy style (already in use in codebase)

**Required Components**:
- One-line summary
- Extended description (for complex functions)
- Parameters section with types
- Returns section with types
- Examples section with doctestable code
- Notes section (optional, for caveats)
- References section (optional, for papers)

**Example**:
```python
def create_signed_lattice(
    side: int,
    pflip: float,
    seed: int = 42
) -> nx.Graph:
    """
    Create a square lattice with randomly flipped edge signs.

    Parameters
    ----------
    side : int
        Number of nodes along each side of the square lattice.
    pflip : float
        Probability of flipping an edge from positive to negative.
    seed : int, optional
        Random seed for reproducibility. Default is 42.

    Returns
    -------
    nx.Graph
        A NetworkX graph representing the signed lattice with 'sign'
        edge attribute (+1 or -1).

    Examples
    --------
    >>> import networkx as nx
    >>> G = create_signed_lattice(side=4, pflip=0.2, seed=12345)
    >>> G.number_of_nodes()
    16
    >>> G.number_of_edges()
    24
    >>> all('sign' in G.edges[e] for e in G.edges)
    True

    Notes
    -----
    This function uses NetworkX's grid_2d_graph internally and then
    randomly assigns negative signs to edges based on pflip.

    See Also
    --------
    nx_patches.Lattice2D : Full-featured lattice class
    """
    # Implementation...
```

### Code Example Standards

**Requirements**:
1. All examples must be reproducible with fixed seeds
2. Examples should be minimal (focus on one concept)
3. Examples must pass `make doctest`
4. Use `>>>` for interactive examples
5. Show output when informative

**Doctest Directives**:
- Use `# doctest: +SKIP` for optional dependencies (cupy, graph-tool)
- Use `# doctest: +ELLIPSIS` for truncated output
- Use `# doctest: +NORMALIZE_WHITESPACE` when needed

### Testing Policy

**Requirements**:
1. All code examples in documentation must pass doctest
2. Examples must use deterministic operations or fixed seeds
3. C programs should have minimal Python wrapper examples
4. Long-running examples should be marked as `+SKIP` with note

---

## Dependencies and Requirements

### New Documentation Dependencies

Add to `pyproject.toml` under `[project.optional-dependencies]`:

```toml
[project.optional-dependencies]
docs = [
    "sphinx>=7.0.0",
    "sphinx-rtd-theme>=2.0.0",
    "sphinx-autodoc-typehints>=1.24.0",
    "myst-parser>=2.0.0",
    "nbsphinx>=0.9.0",
    "sphinx-multiversion>=0.2.4",
    "sphinx-copybutton>=0.5.0",
    "sphinxcontrib-bibtex>=2.6.0",  # For references
]
```

### Installation Command

```bash
pip install -e ".[docs]"
```

---

## Build and Test Commands

### Documentation Build Commands

```bash
# Build HTML documentation
make -C docs html

# Build PDF documentation
make -C docs latexpdf

# Run doctest validation
make -C docs doctest

# Live reload server for development
make -C docs livehtml

# Clean build artifacts
make -C docs clean

# Build all formats
make -C docs all
```

### Integration with Main Test Suite

Add to main project workflow:

```bash
# Run all tests including doctests
pytest --doctest-modules src/

# Run only doctests
pytest --doctest-modules --ignore=test/ src/

# Quick validation
make -C docs doctest
```

---

## Continuous Integration

### GitHub Actions Workflow

**Workflow**: `.github/workflows/docs.yml`

**Triggers**:
- Push to main branch
- Pull requests
- Manual workflow dispatch
- Release tags

**Jobs**:
1. Build documentation
2. Run doctest validation
3. Deploy to GitHub Pages (on main branch only)

**Artifacts**:
- HTML documentation
- Doctest results
- Coverage report (docstring coverage)

---

## Success Metrics

### Quantitative Goals

- [ ] 100% of public functions have docstrings
- [ ] 90%+ of public functions have examples
- [ ] 100% of documentation examples pass doctest
- [ ] Zero Sphinx build warnings
- [ ] < 5 minute documentation build time
- [ ] Documentation accessible at stable URL

### Qualitative Goals

- [ ] New users can install and run first example in < 10 minutes
- [ ] Developers can understand C/C++ integration from docs alone
- [ ] Documentation is searchable and navigable
- [ ] Examples are copy-pasteable and work immediately
- [ ] Theory section connects to academic literature

---

## Timeline and Milestones

**Phase 1: Setup** (1-2 sessions)
- Milestone: Documentation builds and deploys automatically

**Phase 2: API Reference** (2-3 sessions)
- Milestone: Complete API reference with validated examples

**Phase 3: User Guide** (2-3 sessions)
- Milestone: Comprehensive tutorials and examples

**Phase 4: Developer Docs** (1-2 sessions)
- Milestone: Complete developer guide

**Phase 5: Polish** (1 session)
- Milestone: Production-ready documentation

**Total Estimated Effort**: 7-11 working sessions

---

## Risks and Mitigation

### Risk: Existing code lacks docstrings
**Mitigation**: Prioritize most-used modules first, add docstrings incrementally

### Risk: C/C++ code is hard to document
**Mitigation**: Focus on Python wrappers, document interfaces not internals

### Risk: Examples may break with dependency updates
**Mitigation**: Pin documentation build dependencies, run doctest in CI

### Risk: Large scope may delay completion
**Mitigation**: Phased approach, each phase delivers value independently

---

## Open Questions

1. **Hosting preference**: GitHub Pages or ReadTheDocs?
   - **Recommendation**: GitHub Pages (simpler, already using GitHub)

2. **Version strategy**: Document all versions or only latest + stable?
   - **Recommendation**: Start with latest only, add versioning in Phase 5

3. **Notebook inclusion**: Which notebooks to include in docs?
   - **Recommendation**: Create tutorial-specific notebooks, reference research notebooks

4. **Theory depth**: How much mathematical detail?
   - **Recommendation**: Overview with references to papers, not full derivations

5. **C/C++ API**: Document C functions directly or only Python wrappers?
   - **Recommendation**: Python wrappers primary, C internals in developer guide

---

## Next Steps

1. **User approval of this plan**
2. **Begin Phase 1: Setup and Infrastructure**
3. **Iterate with user feedback after each phase**

---

## References

- Sphinx documentation: https://www.sphinx-doc.org/
- NumPy docstring guide: https://numpydoc.readthedocs.io/
- ReadTheDocs tutorial: https://docs.readthedocs.io/
- Example Python projects with excellent docs:
  - NumPy: https://numpy.org/doc/
  - NetworkX: https://networkx.org/documentation/
  - SciPy: https://docs.scipy.org/
