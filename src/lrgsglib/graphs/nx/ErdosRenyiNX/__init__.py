"""
ErdosRenyiNX: NetworkX-based Erdos-Renyi signed graph.

This module provides the ErdosRenyiNX class for generating signed
Erdos-Renyi random graphs, extracted to the largest connected component.
"""

from .ErdosRenyiNX import ErdosRenyiNX

# Backward compatibility alias
ErdosRenyi = ErdosRenyiNX

__all__ = [
    "ErdosRenyiNX",
    "ErdosRenyi",  # alias for backward compatibility
]
