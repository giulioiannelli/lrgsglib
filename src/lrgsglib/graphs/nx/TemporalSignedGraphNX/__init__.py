"""
TemporalSignedGraphNX: Deprecated - use TemporalGraphNX directly.

TemporalGraphNX now inherits from SignedGraphNX and includes all
sign-flipping functionality that was previously in this class.
"""

from ..TemporalGraphNX.TemporalGraphNX import TemporalGraphNX

# Alias for backward compatibility (no deprecation warning on import,
# only the .py module-level warning fires if imported directly)
TemporalSignedGraphNX = TemporalGraphNX

# Backward compatibility alias
TemporalSignedGraph = TemporalGraphNX

__all__ = [
    "TemporalSignedGraphNX",
    "TemporalSignedGraph",  # alias for backward compatibility
]
