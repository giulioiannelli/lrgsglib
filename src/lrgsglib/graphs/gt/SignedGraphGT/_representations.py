"""Representation stubs for SignedGraphGT (NX multi-repr compat).

NX has a multi-representation system (gr['G'], gr['H']) with methods
to update matrices, edge sets, and node/edge mappings.  GT has a single
representation.  These stubs ensure engine-agnostic code that calls
upd_* or zip_* does not crash.
"""

from __future__ import annotations

from typing import Any, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from .SignedGraphGT import SignedGraphGT


def upd_graph_matrices(
    self: "SignedGraphGT", on_g: Optional[str] = None, **kwargs: Any
) -> None:
    """Invalidate cached matrices (GT equivalent of NX upd_graph_matrices)."""
    self.invalidate_cache()


def upd_edge_sets(
    self: "SignedGraphGT", on_g: Optional[str] = None, **kwargs: Any
) -> None:
    """No-op (GT tracks signs in property maps, not separate sets)."""
    pass


def upd_Degree(
    self: "SignedGraphGT", on_g: Optional[str] = None, **kwargs: Any
) -> None:
    """No-op (GT computes degree on the fly)."""
    pass


def upd_GraphRepr_All(
    self: "SignedGraphGT",
    on_g: Optional[str] = None,
    also_itself: bool = True,
    **kwargs: Any,
) -> None:
    """Invalidate cache (GT has single representation)."""
    self.invalidate_cache()


def upd_NodeMap(
    self: "SignedGraphGT", on_g: Optional[str] = None, **kwargs: Any
) -> None:
    """No-op (GT uses integer vertex indices directly)."""
    pass


def upd_EdgeMap(
    self: "SignedGraphGT", on_g: Optional[str] = None, **kwargs: Any
) -> None:
    """No-op (GT edges are graph-tool Edge objects)."""
    pass


def upd_ReprMaps(
    self: "SignedGraphGT", on_g: Optional[str] = None, **kwargs: Any
) -> None:
    """No-op (GT has no representation mappings)."""
    pass


def zip_reprNodes(
    self: "SignedGraphGT", x: Any, on_g: Optional[str] = None
) -> dict:
    """Identity mapping (GT: node indices are their own map)."""
    if hasattr(x, "__iter__"):
        return {i: v for i, v in enumerate(x)}
    return {}


def zip_reprEdges(
    self: "SignedGraphGT", x: Any, on_g: Optional[str] = None
) -> dict:
    """Identity mapping (GT: edges are (src, tgt) tuples)."""
    if hasattr(x, "__iter__"):
        return {i: v for i, v in enumerate(x)}
    return {}
