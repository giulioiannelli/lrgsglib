"""
TemporalSignedGraphNX: Temporal network with signed edges that can flip over time.

This module provides the TemporalSignedGraphNX class which extends
TemporalGraphNX with specific support for signed edge dynamics,
including frustration evolution and sign flipping events.
"""

from .TemporalSignedGraphNX import TemporalSignedGraphNX

# Backward compatibility alias
TemporalSignedGraph = TemporalSignedGraphNX

__all__ = [
    "TemporalSignedGraphNX",
    "TemporalSignedGraph",  # alias for backward compatibility
]
