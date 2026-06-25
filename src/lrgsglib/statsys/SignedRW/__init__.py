"""Signed random-walker dynamics family.

Public API
----------
* :class:`SignedWalker` — abstract base for position-tracked walkers.
* :class:`AbsorbingWalker`, :class:`KillingWalker`, :class:`StickyWalker`
  — concrete rule classes (D3 / D2 / D1 in the notebook convention).
* :class:`SignedSpinCopy` — legacy spin-copy-with-sign dynamics
  (previously exported as ``SignedRW``; that name is kept as a
  backwards-compatible alias).
* :mod:`_overlap` — offline overlap kernels between visit-count fields.
"""

from .SignedSpinCopy import SignedSpinCopy
from .SignedWalker import (
    AbsorbingWalker,
    KillingWalker,
    SignedWalker,
    StickyWalker,
)
from . import _overlap
from . import _solvers  # noqa: F401  (registers the SignedRW solvers at import)

# Backward-compatible alias: this model was formerly named `SignedRW`.
SignedRW = SignedSpinCopy

__all__ = [
    "SignedWalker",
    "AbsorbingWalker",
    "KillingWalker",
    "StickyWalker",
    "SignedSpinCopy",
    "SignedRW",  # legacy alias
]
