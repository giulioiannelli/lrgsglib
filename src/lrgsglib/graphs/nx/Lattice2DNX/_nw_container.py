"""Negative-link (``nwDict``) pattern container for :class:`Lattice2DNX`.

A thin specialization of the engine-neutral :class:`NwContainer` that enables the
full pattern set — ``single`` / ``singleXERR`` / ``singleZERR`` / ``rand`` /
``randXERR`` / ``randZERR`` — exactly mirroring :class:`Lattice2DGTnwContainer`.
ZERR cells come from the shared canonical / coordinate-oriented ``cell_edges``
seam (``SignedGraph.cell_edges``), so the NX and graph-tool lattices build the
same cells node-for-node. This replaces the old bespoke per-geometry container
(``get_links_square`` / ``_triangle`` / ``_hexagon``), which broke cross-engine
reproducibility and raised ``KeyError`` on the ``oct_sqr`` / ``kgm`` geometries.
"""
from __future__ import annotations

from ..._shared._nw_container import NwContainer


class Lattice2DNXnwContainer(NwContainer):
    """Full nwDict (central single + random, both XERR and ZERR)."""

    build_single = True
    build_zerr = True


__all__ = ["Lattice2DNXnwContainer"]
