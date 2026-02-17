"""Partitioning / subgraph methods for SignedGraphGT (NX-compatible API).

Methods here supplement the clustering already on SignedGraphGT
(make_clustersYN, get_eigV_cluster_sizes, get_cluster_distribution,
compute_pinf, handle_no_clust) with subgraph extraction, FM/AFM
region detection, and edge-based connected components.
"""

from __future__ import annotations

import logging
from typing import Any, Optional, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components as sp_connected_components

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from .SignedGraphGT import SignedGraphGT

try:
    import graph_tool as gt
    from graph_tool import GraphView, Graph
    from graph_tool.topology import label_components
except ImportError:
    pass


def get_subgraph_from_nodes(
    self: "SignedGraphGT",
    list_of_nodes: list,
    on_g: Optional[str] = None,
) -> "GraphView":
    """Extract a subgraph induced by the given nodes.

    Parameters
    ----------
    list_of_nodes : list[int]
        Vertex indices to include.
    on_g : str, optional
        Ignored (GT single representation).

    Returns
    -------
    GraphView
        Filtered view of the graph containing only the specified nodes.
    """
    vfilt = self.G.new_vertex_property("bool", val=False)
    for n in list_of_nodes:
        vfilt[self.G.vertex(n)] = True
    return GraphView(self.G, vfilt=vfilt)


def get_nodes_subgraph_by_kv(
    self: "SignedGraphGT",
    k: str,
    val: Any,
    on_g: Optional[str] = None,
):
    """Partition nodes by vertex property value and return two subgraphs.

    Parameters
    ----------
    k : str
        Vertex property name to partition on.
    val : number or callable
        If number, nodes where ``prop[v] == val`` go to G_yes.
        If callable, ``val(prop[v])`` determines membership.
    on_g : str, optional
        Ignored (GT single representation).

    Returns
    -------
    tuple[GraphView, GraphView]
        ``(G_yes, G_no)`` — subgraph views for matching / non-matching nodes.
    """
    prop = self.G.vertex_properties.get(k)
    if prop is None:
        raise KeyError(f"No vertex property '{k}' on graph")

    predicate = val if callable(val) else (lambda x, _v=val: x == _v)

    mask_yes = self.G.new_vertex_property("bool", val=False)
    mask_no = self.G.new_vertex_property("bool", val=False)

    for v in self.G.vertices():
        if predicate(prop[v]):
            mask_yes[v] = True
        else:
            mask_no[v] = True

    return GraphView(self.G, vfilt=mask_yes), GraphView(self.G, vfilt=mask_no)


def make_graphYN(
    self: "SignedGraphGT",
    k: str,
    val: Any,
    on_g: Optional[str] = None,
) -> None:
    """Create Y/N subgraph partitions from a vertex property.

    Parameters
    ----------
    k : str
        Vertex property name.
    val : number or callable
        Partition criterion.
    on_g : str, optional
        Ignored.
    """
    g_yes, g_no = get_nodes_subgraph_by_kv(self, k, val, on_g)
    # Store in a simple dict structure for retrieval
    if not hasattr(self, "_graph_yn"):
        self._graph_yn = {}
    key = val if not callable(val) else str(val)
    self._graph_yn[(k, key)] = (g_yes, g_no)


def make_connected_component_by_edge(
    self: "SignedGraphGT",
    edge_attr: str = "sign",
    value: int = 1,
    on_g: Optional[str] = None,
) -> None:
    """Extract largest connected component of edges with a given sign.

    Parameters
    ----------
    edge_attr : str
        Edge property name (default ``'sign'``).
    value : int
        Edge value to filter on (default ``1`` = positive).
    on_g : str, optional
        Ignored.

    Sets
    ----
    self.largest_cc : set[int]
        Node indices in the largest connected component.
    self.largest_cc_subgraph : GraphView
        Filtered graph view of that component.
    """
    prop = self.G.edge_properties.get(edge_attr)
    if prop is None:
        raise KeyError(f"No edge property '{edge_attr}'")

    efilt = self.G.new_edge_property("bool")
    for e in self.G.edges():
        efilt[e] = (prop[e] == value)

    gv = GraphView(self.G, efilt=efilt)

    # Find connected components in the edge-filtered view
    comp, hist = label_components(gv)

    if len(hist) == 0:
        raise ValueError(f"No connected components found with {edge_attr}={value}.")

    giant_label = int(np.argmax(hist))

    # Build vertex filter for the giant component
    vfilt = self.G.new_vertex_property("bool", val=False)
    for v in self.G.vertices():
        if comp[v] == giant_label:
            vfilt[v] = True

    self.largest_cc = {int(v) for v in self.G.vertices() if vfilt[v]}
    self.largest_cc_subgraph = GraphView(self.G, vfilt=vfilt, efilt=efilt)


def get_ferroAntiferro_regions(
    self: "SignedGraphGT",
    attr_str: str = "s",
    on_g: Optional[str] = None,
) -> tuple[list[int], list[int]]:
    """Identify ferromagnetic and antiferromagnetic ordering regions.

    A node is *ferro* if all its neighbours have the same attribute value.
    A node is *anti-ferro* if all neighbours have the opposite value.

    Parameters
    ----------
    attr_str : str
        Vertex property name holding spin values (default ``'s'``).
    on_g : str, optional
        Ignored.

    Returns
    -------
    tuple[list[int], list[int]]
        ``(ferroGroup, antiGroup)``
    """
    prop = self.G.vertex_properties.get(attr_str)
    if prop is None:
        raise KeyError(f"No vertex property '{attr_str}' on graph")

    ferroGroup: list[int] = []
    antiGroup: list[int] = []

    for v in self.G.vertices():
        node_val = prop[v]
        neighbors = list(v.out_neighbors())
        if not neighbors:
            continue
        if all(prop[n] == node_val for n in neighbors):
            ferroGroup.append(int(v))
        elif all(prop[n] == -node_val for n in neighbors):
            antiGroup.append(int(v))

    return ferroGroup, antiGroup
