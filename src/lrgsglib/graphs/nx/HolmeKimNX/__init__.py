"""
HolmeKimNX: NetworkX-based Holme-Kim model (BA with clustering).

This module provides the HolmeKimNX class for generating signed
scale-free networks with tunable clustering.
"""

from .HolmeKimNX import HolmeKimNX

# Backward compatibility alias
HolmeKim = HolmeKimNX

__all__ = [
    "HolmeKimNX",
    "HolmeKim",  # alias for backward compatibility
]
