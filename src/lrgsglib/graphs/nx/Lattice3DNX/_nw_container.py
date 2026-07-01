"""Negative-link (``nwDict``) pattern container for :class:`Lattice3DNX`.

A thin specialization of the engine-neutral :class:`NwContainer` that enables the
central-defect patterns — ``single`` / ``singleXERR`` / ``rand`` / ``randXERR`` —
exactly mirroring :class:`Lattice3DGTnwContainer`. No ZERR loop patterns are built
by default (a 3D lattice node sits on several equal-length faces with no oriented
picker to disambiguate them); a structured ``*ZERR`` support still builds on demand
via the shared ``build_zerr`` override, resolving through the canonical,
engine-independent ``cell_edges`` seam (``SignedGraph.cell_edges``). This replaces
the old bespoke nested container, whose ``random`` reference and hand-rolled XERR
loop duplicated the shared machinery.
"""
from __future__ import annotations

from ..._shared._nw_container import NwContainer


class Lattice3DNXnwContainer(NwContainer):
    """3D-lattice nwDict (central single + random, XERR only)."""

    build_single = True
    build_zerr = False


__all__ = ["Lattice3DNXnwContainer"]
