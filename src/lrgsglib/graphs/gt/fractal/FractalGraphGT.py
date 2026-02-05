"""
FractalGraphGT - Falls back to NX implementation (no native GT implementation).
"""

import warnings
from ...nx.fractal import FractalGraphNX

warnings.warn(
    "FractalGraphGT: No native graph-tool implementation. Using NetworkX backend.",
    stacklevel=2,
)

FractalGraphGT = FractalGraphNX

__all__ = ["FractalGraphGT"]
