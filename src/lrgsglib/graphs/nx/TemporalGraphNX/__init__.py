"""
TemporalGraphNX: Time-varying network with snapshot-based representation.

This module provides the TemporalGraphNX class for analyzing networks that
evolve over time, with support for both discrete snapshots and continuous-time
edge streams.
"""

from .TemporalGraphNX import TemporalGraphNX

# Backward compatibility alias
TemporalGraph = TemporalGraphNX

__all__ = [
    "TemporalGraphNX",
    "TemporalGraph",  # alias for backward compatibility
]
