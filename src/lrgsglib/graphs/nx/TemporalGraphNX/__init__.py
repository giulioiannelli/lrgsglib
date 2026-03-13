"""
TemporalGraphNX: Time-varying signed network inheriting from SignedGraphNX.

This module provides the TemporalGraphNX class for analyzing networks that
evolve over time, with full signed graph functionality (spectral analysis,
entropy, clustering) via inheritance from SignedGraphNX.
"""

from .TemporalGraphNX import TemporalGraphNX

# Backward compatibility alias
TemporalGraph = TemporalGraphNX

__all__ = [
    "TemporalGraphNX",
    "TemporalGraph",  # alias for backward compatibility
]
