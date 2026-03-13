"""Temporal graph facade module.

Re-exports temporal graph classes. TemporalSignedGraphNX is deprecated —
TemporalGraphNX now inherits from SignedGraphNX directly.
"""

from .TemporalGraphNX import TemporalGraphNX

# Backward compatibility aliases
TemporalGraph = TemporalGraphNX
TemporalSignedGraphNX = TemporalGraphNX  # deprecated alias
TemporalSignedGraph = TemporalGraphNX  # deprecated alias

__all__ = [
    "TemporalGraphNX",
    "TemporalSignedGraphNX",
    "TemporalGraph",
    "TemporalSignedGraph",
]
