"""Fractal graph facade module.

Re-exports fractal graph classes from FractalGraphNX, DGMgraphNX, and SierpinskiNX
packages for cleaner imports.
"""

from .DGMgraphNX import (  # Backward compatibility alias
    DGMgraph,
    DGMgraphNX,
)
from .FractalGraphNX import (  # Backward compatibility alias
    FractalGraph,
    FractalGraphNX,
    FractalNwContainerBase,
)
from .SierpinskiNX import (  # Backward compatibility alias
    SierpinskiGraph,
    SierpinskiNX,
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
