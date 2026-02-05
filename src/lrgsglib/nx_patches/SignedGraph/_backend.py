"""Backend module shim for backward compatibility.

DEPRECATED: This module re-exports from lrgsglib.graphs.nx.SignedGraphNX._backend.
Import directly from lrgsglib.graphs.nx.SignedGraphNX._backend instead.
"""

import warnings

warnings.warn(
    "Importing from lrgsglib.nx_patches.SignedGraph._backend is deprecated. "
    "Use lrgsglib.graphs.nx.SignedGraphNX._backend instead.",
    DeprecationWarning,
    stacklevel=2
)

# Re-export everything from the new location
from ...graphs.nx.SignedGraphNX._backend import (
    Backend,
    ArrayBackend,
    BackendManager,
    NumpyBackend,
    ScipyBackend,
    CupyBackend,
)

__all__ = [
    "Backend",
    "ArrayBackend",
    "BackendManager",
    "NumpyBackend",
    "ScipyBackend",
    "CupyBackend",
]
