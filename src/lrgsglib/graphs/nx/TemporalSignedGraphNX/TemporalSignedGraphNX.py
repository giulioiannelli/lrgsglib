"""
TemporalSignedGraphNX: Deprecated - use TemporalGraphNX directly.

TemporalGraphNX now inherits from SignedGraphNX and includes all
sign-flipping functionality that was previously in this class.

This module provides a backward-compatible alias.
"""

from ..TemporalGraphNX.TemporalGraphNX import TemporalGraphNX

TemporalSignedGraphNX = TemporalGraphNX
