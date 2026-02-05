"""Temporal graph facade module.

Re-exports temporal graph classes from TemporalGraphNX and TemporalSignedGraphNX
packages for cleaner imports.
"""

from .TemporalGraphNX import TemporalGraphNX
from .TemporalSignedGraphNX import TemporalSignedGraphNX

# Backward compatibility aliases
TemporalGraph = TemporalGraphNX
TemporalSignedGraph = TemporalSignedGraphNX

__all__ = [
    "TemporalGraphNX",
    "TemporalSignedGraphNX",
    "TemporalGraph",
    "TemporalSignedGraph",
]
