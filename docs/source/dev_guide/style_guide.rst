Style Guide
===========

This document describes code style conventions, documentation requirements,
and best practices for lrgsglib.

Python Style
------------

General Conventions
~~~~~~~~~~~~~~~~~~~

- **PEP 8** - Follow Python style guide
- **Line length** - 80 characters maximum
- **Python version** - 3.12+ features allowed

Formatting Tools
~~~~~~~~~~~~~~~~

lrgsglib uses these formatting tools:

.. code-block:: bash

   # Format code
   black src/ --line-length 80

   # Sort imports
   isort src/ --profile black

   # Check style
   flake8 src/

   # Type checking
   mypy src/

Configuration in ``pyproject.toml``:

.. code-block:: toml

   [tool.black]
   line-length = 80

   [tool.isort]
   profile = "black"
   line_length = 80

Naming Conventions
~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1

   * - Type
     - Convention
     - Example
   * - Functions
     - snake_case
     - ``compute_laplacian_spectrum``
   * - Variables
     - snake_case
     - ``eigenvalues``
   * - Classes
     - PascalCase
     - ``SignedGraph``
   * - Constants
     - UPPER_CASE
     - ``DEFAULT_BACKEND``
   * - Private
     - _leading_underscore
     - ``_internal_helper``
   * - Module-level private
     - _leading_underscore
     - ``_VALID_BACKENDS``

**Descriptive names:**

.. code-block:: python

   # Good
   def prepare_lattice(args):
       ...

   # Bad
   def prep_lat(a):
       ...

Type Hints
----------

All new code must have type hints:

.. code-block:: python

   from typing import Any, Callable, Dict, List, Optional, Tuple, Union
   from numpy.typing import NDArray

   def compute_entropy(
       eigenvalues: NDArray,
       num_nodes: int,
       time_grid: Optional[NDArray] = None,
   ) -> Tuple[NDArray, NDArray]:
       """Compute entropy from eigenvalues."""
       ...

Common patterns:

.. code-block:: python

   # Optional parameter
   def func(value: Optional[str] = None) -> None:
       ...

   # Union types
   def func(value: Union[int, float]) -> float:
       ...

   # Callable
   def func(callback: Callable[[int], NDArray]) -> None:
       ...

   # Self reference (for methods)
   def method(self: "SignedGraph") -> "SignedGraph":
       ...

Docstrings
----------

Use NumPy docstring format:

.. code-block:: python

   def compute_laplacian_spectrum(
       self: "SignedGraph",
       backend: Optional[str] = None,
       keep_sparse: bool = None,
   ) -> None:
       """
       Compute eigenvalues of the signed Laplacian.

       This function computes the eigenvalues of the signed Laplacian matrix
       without computing eigenvectors, which is faster when eigenvectors
       are not needed.

       Parameters
       ----------
       backend : str, optional
           Backend for computation ('numpy', 'scipy', 'cupy').
           Defaults to instance backend.
       keep_sparse : bool, optional
           If True, use sparse methods (computes N-2 eigenvalues).
           If None, automatically decide based on graph size.

       Returns
       -------
       None
           Results stored in instance attribute:

           - self.eigv : ndarray of shape (N,)
               Eigenvalues in ascending order.

       Notes
       -----
       This function is useful when you only need spectral properties
       without the eigenvectors themselves.

       **Memory Usage:**
       Only stores N eigenvalues (8N bytes for float64).

       See Also
       --------
       compute_laplacian_spectrum_weigV : Compute eigenvalues and eigenvectors
       compute_k_eigvV : Compute k smallest eigenvalues

       Examples
       --------
       >>> sg = SignedGraph(G, pflip=0.3, seed=42)
       >>> sg.compute_laplacian_spectrum()
       >>> print(sg.eigv.shape)
       (100,)
       """

Required Sections
~~~~~~~~~~~~~~~~~

1. **One-line summary** - What does it do?
2. **Extended description** - Details (optional)
3. **Parameters** - Input arguments
4. **Returns** - Output values
5. **Notes** - Important details (optional)
6. **See Also** - Related functions (optional)
7. **Examples** - Runnable code (recommended)

Module Organization
-------------------

File Structure
~~~~~~~~~~~~~~

.. code-block:: python

   """
   Module docstring describing the module's purpose.

   This module provides utilities for spectral analysis of signed graphs.
   """

   from __future__ import annotations

   # Standard library
   import functools
   from pathlib import Path
   from typing import Optional

   # Third-party
   import numpy as np
   import scipy.linalg as la

   # Local
   from ..config import LRGSG_DATA
   from ..utils.basic import normalize_array

   __all__ = [
       "public_function",
       "PublicClass",
   ]

   # Module-level constants
   _VALID_BACKENDS = {"numpy", "scipy", "cupy"}
   DEFAULT_TOLERANCE = 1e-10


   def public_function(arg: int) -> float:
       """Public function with docstring."""
       return _private_helper(arg)


   def _private_helper(arg: int) -> float:
       """Private helper function."""
       return float(arg)


   class PublicClass:
       """Public class with docstring."""

       def __init__(self, value: int) -> None:
           self.value = value

Export Control
~~~~~~~~~~~~~~

Always define ``__all__`` to document public API:

.. code-block:: python

   __all__ = [
       "compute_entropy",
       "compute_specific_heat",
       "EntropyResult",
   ]

Function Design
---------------

Single Responsibility
~~~~~~~~~~~~~~~~~~~~~

Each function should do one thing:

.. code-block:: python

   # Good - separate concerns
   def load_data(path: Path) -> NDArray:
       """Load data from file."""
       ...

   def process_data(data: NDArray) -> NDArray:
       """Process loaded data."""
       ...

   def save_results(results: NDArray, path: Path) -> None:
       """Save results to file."""
       ...

   # Bad - mixed concerns
   def load_process_save(input_path: Path, output_path: Path) -> None:
       """Load, process, and save data."""
       ...

Function Length
~~~~~~~~~~~~~~~

Keep functions under 50 lines. Split larger functions:

.. code-block:: python

   # Good - composed from smaller functions
   def analyze_graph(sg: SignedGraph) -> AnalysisResult:
       """Perform complete graph analysis."""
       spectrum = compute_spectrum(sg)
       entropy = compute_entropy(spectrum)
       clusters = find_clusters(sg, spectrum)
       return AnalysisResult(spectrum, entropy, clusters)

Error Handling
~~~~~~~~~~~~~~

Provide helpful error messages:

.. code-block:: python

   def resolve_backend(backend: str) -> str:
       """Validate and normalize backend identifier."""
       backend_norm = backend.lower()
       if backend_norm not in _VALID_BACKENDS:
           raise ValueError(
               f"Unsupported backend '{backend}'. "
               f"Expected one of {_VALID_BACKENDS}."
           )
       if backend_norm == "cupy":
           try:
               import cupy
           except ImportError as exc:
               raise RuntimeError(
                   "Backend 'cupy' requested but CuPy is not available. "
                   "Install cupy or choose 'numpy'."
               ) from exc
       return backend_norm

Class Design
------------

Private Method Pattern
~~~~~~~~~~~~~~~~~~~~~~

For large classes, split methods into separate files:

.. code-block:: text

   graphs/nx/SignedGraphNX/
   ├── SignedGraphNX.py   ← Main class
   ├── _spectral.py       ← Spectral methods
   ├── _dynamics.py       ← Dynamics methods
   └── _loaders.py        ← I/O methods

In main class:

.. code-block:: python

   class SignedGraph:
       # Import methods from private modules
       from ._spectral import (
           compute_laplacian_spectrum,
           compute_k_eigvV,
       )
       from ._dynamics import (
           init_ising_dynamics,
           run_dynamics,
       )

Initialization
~~~~~~~~~~~~~~

.. code-block:: python

   class SignedGraph:
       """Signed graph with positive and negative edges."""

       def __init__(
           self,
           graph: nx.Graph,
           *,  # Force keyword arguments
           pflip: float = 0.0,
           seed: Optional[int] = None,
       ) -> None:
           """
           Initialize signed graph.

           Parameters
           ----------
           graph : nx.Graph
               Base NetworkX graph.
           pflip : float, default 0.0
               Fraction of edges to flip to negative.
           seed : int, optional
               Random seed for reproducibility.
           """
           self._graph = graph
           self._pflip = pflip
           self._seed = seed
           self._eigv: Optional[NDArray] = None

Comments
--------

Comment Philosophy
~~~~~~~~~~~~~~~~~~

- Explain **why**, not **what**
- Code should be self-documenting
- Use comments for non-obvious logic

.. code-block:: python

   # Good - explains why
   # Use sparse solver for large graphs to avoid memory issues
   if N > 10000:
       use_sparse = True

   # Bad - explains what (obvious from code)
   # Set use_sparse to True
   use_sparse = True

TODO Comments
~~~~~~~~~~~~~

.. code-block:: python

   # TODO(username): Fix edge case when graph is disconnected
   # FIXME: This is a temporary workaround for issue #42

Imports
-------

Import Order
~~~~~~~~~~~~

1. Standard library
2. Third-party packages
3. Local imports

.. code-block:: python

   # Standard library
   import os
   from pathlib import Path
   from typing import Optional

   # Third-party
   import numpy as np
   import networkx as nx
   from scipy import linalg

   # Local
   from ..config import LRGSG_DATA
   from ..utils.basic import normalize

Import Style
~~~~~~~~~~~~

.. code-block:: python

   # Good - explicit imports
   from numpy import ndarray
   import scipy.linalg as la

   # Avoid - star imports
   from numpy import *

   # Exception - __init__.py for re-export
   from .module import *

Testing Requirements
--------------------

All new code requires tests:

.. code-block:: python

   def test_compute_entropy_basic():
       """Test entropy computation with known input."""
       eigenvalues = np.array([0.0, 1.0, 2.0, 3.0])
       result = compute_entropy(eigenvalues, num_nodes=4)
       assert result.shape[0] > 0


   @pytest.mark.parametrize("backend", ["numpy", "scipy"])
   def test_compute_entropy_backends(backend):
       """Test entropy computation works with all backends."""
       ...

Performance Considerations
--------------------------

Vectorization
~~~~~~~~~~~~~

Prefer NumPy operations over loops:

.. code-block:: python

   # Good - vectorized
   result = np.exp(-eigenvalues * time)

   # Bad - loop
   result = np.array([np.exp(-e * time) for e in eigenvalues])

Memory Efficiency
~~~~~~~~~~~~~~~~~

.. code-block:: python

   # Good - in-place operations
   data *= 2.0

   # Consider memory for large arrays
   if data.nbytes > 1e9:  # > 1GB
       # Process in chunks
       ...

Avoid Premature Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Write clear code first
- Profile to find bottlenecks
- Optimize only where needed

Anti-Patterns to Avoid
----------------------

1. **Global state** - Use parameters and return values
2. **Magic numbers** - Define constants with names
3. **Deep nesting** - Extract to functions
4. **Long functions** - Split into smaller pieces
5. **Mixed concerns** - Separate I/O from logic
6. **Hardcoded paths** - Use configuration
7. **Missing error handling** - Validate inputs
8. **No type hints** - Always add types
9. **No tests** - Test everything
10. **No docstrings** - Document public API

See Also
--------

- :doc:`architecture` - System design patterns
- :doc:`testing` - Writing tests
- :doc:`contributing` - Development workflow
