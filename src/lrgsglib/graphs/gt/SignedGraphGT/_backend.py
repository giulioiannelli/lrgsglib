"""Backend abstraction layer — re-exports from shared location.

The canonical implementation lives in ``graphs._shared._backend``.
"""

from ..._shared._backend import (  # noqa: F401
    Backend,
    ArrayBackend,
    BackendManager,
    NumpyBackend,
    ScipyBackend,
    CupyBackend,
)
