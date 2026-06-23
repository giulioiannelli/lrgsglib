"""Negative-link (``nwDict``) pattern container for :class:`Lattice3DGT`.

Mirrors ``Lattice3DNX.nwContainer``: a central single defect plus random XERR
stars (``single`` / ``singleXERR`` / ``rand`` / ``randXERR``). No ZERR — the 3D
NX container does not define loop patterns either.
"""
from __future__ import annotations

from ..._shared._nw_container import GTnwContainer


class Lattice3DGTnwContainer(GTnwContainer):
    """3D-lattice nwDict (central single + random, XERR only)."""

    build_single = True
    build_zerr = False


__all__ = ["Lattice3DGTnwContainer"]
