"""
Converters: Bidirectional conversion between NetworkX and graph-tool.

This module provides utilities for converting graphs between NetworkX
and graph-tool formats, preserving edge signs and other attributes.

Functions
---------
nx_to_gt : Convert NetworkX graph to graph-tool Graph
gt_to_nx : Convert graph-tool Graph to NetworkX graph

Example
-------
>>> import networkx as nx
>>> from lrgsglib.graphs.gt._converters import nx_to_gt, gt_to_nx
>>>
>>> # Create NX graph with signs
>>> G_nx = nx.complete_graph(10)
>>> for u, v in G_nx.edges():
...     G_nx[u][v]['sign'] = 1 if (u + v) % 2 == 0 else -1
>>>
>>> # Convert to graph-tool
>>> G_gt = nx_to_gt(G_nx)
>>>
>>> # Convert back
>>> G_nx_back = gt_to_nx(G_gt)
"""

from typing import Any, Optional

import numpy as np

try:
    import graph_tool as gt
    from graph_tool import Graph
    GT_AVAILABLE = True
except ImportError:
    GT_AVAILABLE = False
    Graph = None

import networkx as nx


class GTNXConverter:
    """
    Bidirectional converter between NetworkX and graph-tool graphs.

    This class handles the conversion of graphs including:
    - Node/edge mapping
    - Edge signs (stored as edge property)
    - Other edge/node attributes

    Attributes
    ----------
    node_map_nx_to_gt : dict
        Mapping from NX node IDs to GT vertex indices.
    node_map_gt_to_nx : dict
        Mapping from GT vertex indices to NX node IDs.

    Examples
    --------
    >>> converter = GTNXConverter()
    >>> G_gt = converter.nx_to_gt(G_nx)
    >>> G_nx = converter.gt_to_nx(G_gt)
    """

    def __init__(self):
        self.node_map_nx_to_gt = {}
        self.node_map_gt_to_nx = {}

    def nx_to_gt(
        self,
        G_nx: nx.Graph,
        sign_attr: str = "sign",
        weight_attr: Optional[str] = None,
    ) -> "Graph":
        """
        Convert NetworkX graph to graph-tool Graph.

        Parameters
        ----------
        G_nx : nx.Graph
            Input NetworkX graph.
        sign_attr : str, default "sign"
            Edge attribute name for signs in NX graph.
        weight_attr : str, optional
            Edge attribute name for weights (if any).

        Returns
        -------
        Graph
            graph-tool Graph with edge properties for signs.
        """
        if not GT_AVAILABLE:
            raise ImportError("graph-tool is not installed")

        # Create GT graph
        G_gt = Graph(directed=G_nx.is_directed())

        # Create node mapping
        self.node_map_nx_to_gt = {}
        self.node_map_gt_to_nx = {}

        # Add vertices
        for i, node in enumerate(G_nx.nodes()):
            v = G_gt.add_vertex()
            self.node_map_nx_to_gt[node] = int(v)
            self.node_map_gt_to_nx[int(v)] = node

        # Create edge property maps
        sign_prop = G_gt.new_edge_property("int")
        if weight_attr:
            weight_prop = G_gt.new_edge_property("double")

        # Add edges with properties
        for u, v, data in G_nx.edges(data=True):
            u_gt = self.node_map_nx_to_gt[u]
            v_gt = self.node_map_nx_to_gt[v]
            e = G_gt.add_edge(u_gt, v_gt)

            # Set sign (default +1 if not present)
            sign_prop[e] = data.get(sign_attr, 1)

            if weight_attr and weight_attr in data:
                weight_prop[e] = data[weight_attr]

        # Store property maps
        G_gt.edge_properties["sign"] = sign_prop
        if weight_attr:
            G_gt.edge_properties["weight"] = weight_prop

        # Store node mapping as graph property
        G_gt.graph_properties["node_map_gt_to_nx"] = G_gt.new_graph_property(
            "object", self.node_map_gt_to_nx
        )

        return G_gt

    def gt_to_nx(
        self,
        G_gt: "Graph",
        sign_prop: str = "sign",
        weight_prop: Optional[str] = None,
    ) -> nx.Graph:
        """
        Convert graph-tool Graph to NetworkX graph.

        Parameters
        ----------
        G_gt : Graph
            Input graph-tool Graph.
        sign_prop : str, default "sign"
            Edge property name for signs.
        weight_prop : str, optional
            Edge property name for weights.

        Returns
        -------
        nx.Graph
            NetworkX graph with edge attributes.
        """
        if not GT_AVAILABLE:
            raise ImportError("graph-tool is not installed")

        # Create NX graph
        if G_gt.is_directed():
            G_nx = nx.DiGraph()
        else:
            G_nx = nx.Graph()

        # Try to get stored node mapping
        if "node_map_gt_to_nx" in G_gt.graph_properties:
            node_map = G_gt.graph_properties["node_map_gt_to_nx"]
        else:
            # Default: use vertex indices as node IDs
            node_map = {int(v): int(v) for v in G_gt.vertices()}

        # Add nodes
        G_nx.add_nodes_from(node_map.values())

        # Get edge properties
        has_sign = sign_prop in G_gt.edge_properties
        has_weight = weight_prop and weight_prop in G_gt.edge_properties

        # Add edges
        for e in G_gt.edges():
            u_nx = node_map[int(e.source())]
            v_nx = node_map[int(e.target())]

            edge_data = {}
            if has_sign:
                edge_data["sign"] = G_gt.edge_properties[sign_prop][e]
            if has_weight:
                edge_data["weight"] = G_gt.edge_properties[weight_prop][e]

            G_nx.add_edge(u_nx, v_nx, **edge_data)

        return G_nx


# Module-level convenience functions
_default_converter = None


def _get_converter() -> GTNXConverter:
    """Get or create default converter instance."""
    global _default_converter
    if _default_converter is None:
        _default_converter = GTNXConverter()
    return _default_converter


def nx_to_gt(
    G_nx: nx.Graph,
    sign_attr: str = "sign",
    weight_attr: Optional[str] = None,
) -> "Graph":
    """
    Convert NetworkX graph to graph-tool Graph.

    Parameters
    ----------
    G_nx : nx.Graph
        Input NetworkX graph.
    sign_attr : str, default "sign"
        Edge attribute name for signs.
    weight_attr : str, optional
        Edge attribute name for weights.

    Returns
    -------
    Graph
        graph-tool Graph.

    Examples
    --------
    >>> import networkx as nx
    >>> from lrgsglib.graphs.gt._converters import nx_to_gt
    >>> G_nx = nx.complete_graph(5)
    >>> G_gt = nx_to_gt(G_nx)
    """
    return _get_converter().nx_to_gt(G_nx, sign_attr, weight_attr)


def gt_to_nx(
    G_gt: "Graph",
    sign_prop: str = "sign",
    weight_prop: Optional[str] = None,
) -> nx.Graph:
    """
    Convert graph-tool Graph to NetworkX graph.

    Parameters
    ----------
    G_gt : Graph
        Input graph-tool Graph.
    sign_prop : str, default "sign"
        Edge property name for signs.
    weight_prop : str, optional
        Edge property name for weights.

    Returns
    -------
    nx.Graph
        NetworkX graph.

    Examples
    --------
    >>> from lrgsglib.graphs.gt._converters import gt_to_nx
    >>> G_nx = gt_to_nx(G_gt)
    """
    return _get_converter().gt_to_nx(G_gt, sign_prop, weight_prop)


def get_adjacency_matrix_gt(G_gt: "Graph", signed: bool = True) -> np.ndarray:
    """
    Get adjacency matrix from graph-tool Graph.

    Parameters
    ----------
    G_gt : Graph
        Input graph-tool Graph.
    signed : bool, default True
        If True, use edge signs; otherwise binary adjacency.

    Returns
    -------
    np.ndarray
        Adjacency matrix.
    """
    if not GT_AVAILABLE:
        raise ImportError("graph-tool is not installed")

    from graph_tool.spectral import adjacency

    if signed and "sign" in G_gt.edge_properties:
        # Get signed adjacency
        A = adjacency(G_gt, weight=G_gt.edge_properties["sign"]).toarray()
    else:
        A = adjacency(G_gt).toarray()

    return A


def get_laplacian_matrix_gt(G_gt: "Graph", signed: bool = True) -> np.ndarray:
    """
    Get Laplacian matrix from graph-tool Graph.

    Parameters
    ----------
    G_gt : Graph
        Input graph-tool Graph.
    signed : bool, default True
        If True, use signed Laplacian L = D - A_signed.

    Returns
    -------
    np.ndarray
        Laplacian matrix.
    """
    if not GT_AVAILABLE:
        raise ImportError("graph-tool is not installed")

    from graph_tool.spectral import laplacian

    if signed and "sign" in G_gt.edge_properties:
        # Signed Laplacian: L = D - A_signed
        L = laplacian(G_gt, weight=G_gt.edge_properties["sign"]).toarray()
    else:
        L = laplacian(G_gt).toarray()

    return L
