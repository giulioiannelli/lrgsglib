"""Topology helpers for SignedGraphGT (NX-compatible API).

Methods already directly on SignedGraphGT (get_signed_laplacian,
get_adjacency_matrix, get_degree_matrix, get_abs_degree_matrix,
get_laplacian, get_nodes_list, get_signed_laplacian_embedding,
make_rescaled_signed_laplacian) are NOT duplicated here.
"""

from __future__ import annotations

from typing import Any, Optional, TYPE_CHECKING

import numpy as np

logger = __import__("logging").getLogger(__name__)

if TYPE_CHECKING:
    from .SignedGraphGT import SignedGraphGT


# ------------------------------------------------------------------
# Node / edge accessors
# ------------------------------------------------------------------

def nodes_in(self: "SignedGraphGT", on_g: Optional[str] = None) -> list[int]:
    """Get list of node indices (alias for get_nodes_list)."""
    return self.get_nodes_list()


def get_node_attributes(
    self: "SignedGraphGT",
    attr: str = "pos",
    on_g: Optional[str] = None,
) -> dict:
    """Get vertex property values as a dict.

    Parameters
    ----------
    attr : str
        Vertex property name.
    on_g : str, optional
        Ignored (GT single representation).

    Returns
    -------
    dict
        ``{node_index: value}`` for all vertices.
    """
    prop = self.G.vertex_properties.get(attr)
    if prop is None:
        return {}
    return {int(v): prop[v] for v in self.G.vertices()}


def get_edge_data(
    self: "SignedGraphGT",
    u: int,
    v: int,
    thedata: str = "sign",
    on_g: Optional[str] = None,
) -> Any:
    """Get a specific edge property value.

    Parameters
    ----------
    u, v : int
        Vertex indices.
    thedata : str
        Edge property name (default ``'sign'``).
    on_g : str, optional
        Ignored.

    Returns
    -------
    value
        The edge property value.

    Raises
    ------
    KeyError
        If edge doesn't exist or property not found.
    """
    e = self.G.edge(u, v)
    if e is None:
        raise KeyError(f"Edge ({u}, {v}) not found")
    prop = self.G.edge_properties.get(thedata)
    if prop is None:
        raise KeyError(f"No edge property '{thedata}'")
    return prop[e]


def get_edge_color(
    self: "SignedGraphGT",
    pec: str = "blue",
    nec: str = "red",
    thedata: str = "sign",
    on_g: Optional[str] = None,
    continuous: bool = False,
    cmap: str = "coolwarm",
) -> list:
    """Get edge colors based on sign (or continuous weight).

    Parameters
    ----------
    pec : str
        Positive edge color.
    nec : str
        Negative edge color.
    thedata : str
        Edge property to colour by (default ``'sign'``).
    on_g : str, optional
        Ignored.
    continuous : bool
        If True, normalise values to [0,1] for a colormap.
    cmap : str
        Matplotlib colormap name (only used if ``continuous=True``).

    Returns
    -------
    list
        List of colours, one per edge.
    """
    prop = self.G.edge_properties.get(thedata)
    if prop is None:
        return [pec] * self.num_edges

    values = [prop[e] for e in self.G.edges()]

    if continuous:
        arr = np.array(values, dtype=float)
        vmin, vmax = arr.min(), arr.max()
        if vmax != vmin:
            return list((arr - vmin) / (vmax - vmin))
        return [0.5] * len(arr)

    return [nec if v == -1 else pec for v in values]


def get_graph_neighbors(
    self: "SignedGraphGT",
    node: int,
    on_g: Optional[str] = None,
) -> list[int]:
    """Get neighbours of a node (alias for get_neighbor_indices)."""
    return self.get_neighbor_indices(node)


# ------------------------------------------------------------------
# Multi-representation matrix accessors (GT: single repr, thin wrappers)
# ------------------------------------------------------------------

def get_adjacency_matrix_for(
    self: "SignedGraphGT", on_g: Optional[str] = None
) -> np.ndarray:
    """Get signed adjacency matrix (GT single representation)."""
    return self.get_signed_adjacency()


def get_degree_matrix_for(
    self: "SignedGraphGT", on_g: Optional[str] = None
) -> np.ndarray:
    """Get degree matrix (GT single representation)."""
    return self.get_degree_matrix()


def get_laplacian_matrix_for(
    self: "SignedGraphGT", on_g: Optional[str] = None
) -> np.ndarray:
    """Get Laplacian matrix (GT single representation)."""
    return self.get_laplacian()


def get_signed_degree_matrix_for(
    self: "SignedGraphGT", on_g: Optional[str] = None
) -> np.ndarray:
    """Get absolute degree matrix (GT single representation)."""
    return self.get_abs_degree_matrix()


def get_signed_laplacian_matrix_for(
    self: "SignedGraphGT", on_g: Optional[str] = None
) -> np.ndarray:
    """Get signed Laplacian matrix (GT single representation)."""
    return self.get_signed_laplacian()
