"""Negative-link (``nwDict``) pattern container for :class:`Lattice2DGT`.

A thin specialization of :class:`GTnwContainer` that enables the full pattern
set — ``single`` / ``singleXERR`` / ``singleZERR`` / ``rand`` / ``randXERR`` /
``randZERR`` — mirroring ``Lattice2DNX.nwContainer``. ZERR generalizes to the
lattice's elementary face (3-cycle triangular/kagome/tri_hex, 4-cycle
square/octagon, 6-cycle honeycomb) via the shared geometry helpers.
"""
from __future__ import annotations

from ..._shared._nw_container import GTnwContainer


class Lattice2DGTnwContainer(GTnwContainer):
    """Full nwDict (central single + random, both XERR and ZERR)."""

    build_single = True
    build_zerr = True


__all__ = ["Lattice2DGTnwContainer"]
