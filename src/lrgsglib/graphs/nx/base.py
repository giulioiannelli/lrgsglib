"""Base signed graph facade module.

Re-exports SignedGraphNX class from SignedGraphNX package for cleaner imports.
"""

from .SignedGraphNX import (  # Backward compatibility alias
    ArrayBackend,
    Backend,
    BackendManager,
    SignedGraphNX,
)

# Backward compatibility alias
SignedGraph = SignedGraphNX

__all__ = [
    "SignedGraphNX",
    "SignedGraph",
    "Backend",
    "BackendManager",
    "ArrayBackend",
]
