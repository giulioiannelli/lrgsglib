"""
LFRBenchmarkNX: NetworkX-based LFR benchmark graph.

This module provides the LFRBenchmarkNX class for generating signed
LFR benchmark graphs with planted community structure.
"""

from .LFRBenchmarkNX import LFRBenchmarkNX

# Backward compatibility alias
LFRBenchmark = LFRBenchmarkNX

__all__ = [
    "LFRBenchmarkNX",
    "LFRBenchmark",  # alias for backward compatibility
]
