"""Fractal graph facade module.

Re-exports fractal graph classes from FractalGraphNX, DGMgraphNX, and SierpinskiNX
packages for cleaner imports.
"""

from .FractalGraphNX import (
    FractalGraphNX,
    FractalNwContainerBase,
    # Backward compatibility alias
    FractalGraph,
)
from .DGMgraphNX import (
    DGMgraphNX,
    # Backward compatibility alias
    DGMgraph,
)
from .SierpinskiNX import (
    SierpinskiNX,
    # Backward compatibility alias
    SierpinskiGraph,
)

# Alias for backward compatibility with old naming
SierpinskiGraphNX = SierpinskiNX

__all__ = [
    "FractalGraphNX",
    "FractalGraph",
    "FractalNwContainerBase",
    "DGMgraphNX",
    "DGMgraph",
    "SierpinskiNX",
    "SierpinskiGraphNX",
    "SierpinskiGraph",
]
