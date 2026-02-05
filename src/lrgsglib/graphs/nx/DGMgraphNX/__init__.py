"""
DGMgraphNX: NetworkX-based Dorogovtsev-Goltsev-Mendes signed graph.

This module provides the DGMgraphNX class for generating signed
Dorogovtsev-Goltsev-Mendes graphs, a deterministic scale-free network.
"""

from .DGMgraphNX import DGMgraphNX

# Backward compatibility alias
DGMgraph = DGMgraphNX

__all__ = [
    "DGMgraphNX",
    "DGMgraph",  # alias for backward compatibility
]
