"""Backend abstraction layer — re-exports from shared location.

The canonical implementation lives in ``graphs._shared._backend``.
This module re-exports for backward compatibility with existing NX imports.
"""

from ..._shared._backend import (  # noqa: F401
    ArrayBackend,
    Backend,
    BackendManager,
    CupyBackend,
    NumpyBackend,
    ScipyBackend,
)

__all__ = [
    "Backend",
    "ArrayBackend",
    "BackendManager",
    "NumpyBackend",
    "ScipyBackend",
    "CupyBackend",
]
