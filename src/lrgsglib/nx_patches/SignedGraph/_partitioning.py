import logging
import warnings
from typing import Optional, TYPE_CHECKING

import numpy as np
import networkx as nx

from ...config.const import SG_REPR, SG_WARNMSG_NOCLUST
from ...config.errwar import NoClustError, SignedGraphWarning
from ...utils.tools import ConditionalPartitioning

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from .SignedGraph import SignedGraph


def get_subgraph_from_nodes(
        self: "SignedGraph",
        list_of_nodes, 
        on_g: str = SG_REPR
):
    return self.gr[on_g].subgraph(list_of_nodes)


def get_nodes_subgraph_by_kv(
        self: "SignedGraph",
        k, 
        val, 
        on_g: str = SG_REPR
):
    """
    Optimized version: partition nodes first, then create subgraphs.
    This avoids copying the entire graph and removing nodes one by 
    one.
    """
    G = self.gr[on_g]
    predicate = val if callable(val) else lambda x: x == val
    
    # Partition nodes into two sets based on predicate
    nodes_no = []  # nodes where predicate is True
    nodes_yes = []  # nodes where predicate is False
    
    for node, v in G.nodes(data=k):
        if predicate(v):
            nodes_no.append(node)
        else:
            nodes_yes.append(node)
    
    # Create subgraphs using .subgraph() which is much faster
    # Note: subgraph() returns a view, use .copy() if you need 
    # independent graphs
    G_yes = G.subgraph(nodes_yes).copy()
    G_no = G.subgraph(nodes_no).copy()
    
    return G_yes, G_no


def get_eigV_cluster_sizes(
        self: "SignedGraph",
        which: int = 0, 
        binarize: bool = True, 
        val: ConditionalPartitioning = None,
        on_g: str = SG_REPR, 
        backend: str = 'nx'
):
    """
    Get cluster sizes for a given eigenvector.
    
    Parameters
    ----------
    which : int, default 0
        Index of the eigenvector to analyze.
    binarize : bool, default True
        Whether to binarize the eigenvector values.
    val : ConditionalPartitioning, optional
        Conditional partitioning object for clustering. 
        Defaults to +1.
    on_g : str, default SG_REPR
        Graph representation to use.
    backend : str, default 'nx'
        Backend to use for clustering ('nx' for NetworkX). Other backends
        are the same as in make_clustersYN, i.e. 'scipy' or 'cupy'.
        
    Returns
    -------
    list
        Sorted list of cluster sizes in descending order.
    """
    node_has_attr = all(
        f"eigV{which}" in self.gr[on_g].nodes[node] 
        for node in self.gr[on_g].nodes
    )
    if not node_has_attr:
        self.load_eigV_on_graph(which, on_g, binarize)
    
    # Default to +1 if val not provided
    if val is None:
        val = ConditionalPartitioning(
            key="+1", 
            cond_func=lambda x: x == +1
        )
    
    if not hasattr(self, "clustersY"):
        make_clustersYN(self, f"eigV{which}", val, on_g, backend)
    cl_len = sorted(map(len, self.clustersY), reverse=True)
    return cl_len


def get_cluster_distribution(
        self: "SignedGraph",
        which: int = 0, 
        on_g: str = SG_REPR,
        binarize: bool = True,
        val: ConditionalPartitioning = None
):
    """
    Get distribution of cluster sizes for a given eigenvector.
    
    Parameters
    ----------
    which : int, default 0
        Index of the eigenvector to analyze.
    on_g : str, default SG_REPR
        Graph representation to use.
    binarize : bool, default True
        Whether to binarize the eigenvector values.
    val : ConditionalPartitioning, optional
        Conditional partitioning object for clustering. If None, defaults 
        to +1.
        
    Returns
    -------
    dict
        Dictionary mapping cluster sizes to their frequency.
    """
    cl_len = get_eigV_cluster_sizes(
        self, which, binarize, val, on_g
    )
    dictdist_cluster_sizes = {
        size: cl_len.count(size) for size in set(cl_len)
    }
    return dictdist_cluster_sizes


def get_ferroAntiferro_regions(
        self: "SignedGraph",
        attr_str: str = 's', 
        on_g: str = SG_REPR
):
    antiGroup = []
    ferroGroup = []
    graph = self.gr[on_g]
    for node, att in graph.nodes(data=attr_str):
        neighbors = graph.neighbors(node)
        if all(graph.nodes[n][attr_str] == -att for n in neighbors):
            antiGroup.append(node)
        elif all(graph.nodes[n][attr_str] == att for n in neighbors):
            ferroGroup.append(node)
    return ferroGroup, antiGroup


def make_graphYN(
    self: "SignedGraph",
    k,
    val: ConditionalPartitioning,
    on_g: str = SG_REPR
):
    subgraph_result = self.get_nodes_subgraph_by_kv(k, val.cond_func, on_g)
    self.graph_clustering_utility[k][val.key][on_g] = subgraph_result


def make_clustersYN(
    self: "SignedGraph",
    k,
    val: ConditionalPartitioning,
    on_g: str = SG_REPR,
    backend: str = 'nx',
):
    """Populate Y/N clusters, optionally using SciPy or CuPy connected components."""
    try:
        graphY, graphN = self.graph_clustering_utility[k][val.key][on_g]
    except Exception as e:
        logger.error(
            f"""
                Error occurred while accessing graph clustering utility: {e}
                Recomputing Y/N graphs for k={k}, val={val.key}, on_g={on_g}.
            """
        )
        make_graphYN(self, k, val, on_g)
        graphY, graphN = self.graph_clustering_utility[k][val.key][on_g]

    backend = backend.lower()
    use_gpu = backend == 'cupy'
    xp = np

    def _components(G):
        if backend == 'nx':
            return list(nx.connected_components(G))

        nodes = list(G.nodes())
        if not nodes:
            return []

        adj = nx.to_scipy_sparse_array(G, format='csr')
        if adj.nnz == 0:
            return [{node} for node in nodes]

        if backend == 'scipy':
            from scipy.sparse.csgraph import connected_components as cs_cc

            n_comp, labels = cs_cc(adj, directed=False, return_labels=True)
            labels = np.asarray(labels, dtype=np.int64)
        elif backend == 'cupy':
            import cupy as cp
            from cupyx.scipy.sparse import csr_matrix as cx_csr
            from cupyx.scipy.sparse import csgraph as cx_csgraph

            n_comp, labels = cx_csgraph.connected_components(
                cx_csr(adj), directed=False
            )
            xp_local = cp
            labels = xp_local.asarray(labels, dtype=xp_local.int64)
            labels = cp.asnumpy(labels)
        else:
            raise ValueError("backend must be one of {'nx', 'scipy', 'cupy'}")

        clusters = [set() for _ in range(int(n_comp))]
        for node, lab in zip(nodes, labels):
            clusters[int(lab)].add(node)
        return clusters

    self.clustersY = sorted(_components(graphY), key=len, reverse=True)
    self.clustersN = _components(graphN)

    self.numClustersY = len(self.clustersY)
    self.numClustersN = len(self.clustersN)

    if not self.clustersY:
        warnings.warn(SG_WARNMSG_NOCLUST, SignedGraphWarning)
        self.biggestClSet = []
        self.numClustersBig = 0
        self.gc = None
        return

    self.biggestClSet = self.clustersY
    self.numClustersBig = len(self.biggestClSet)
    self.gc = self.clustersY[0]


def make_eigVclustersYN(
    self: "SignedGraph",
    val: ConditionalPartitioning,
    which: int = 0,
    on_g: str = SG_REPR,
    binarize: bool = True,
    backend: str = 'nx',
):
    self.load_eigV_on_graph(which, on_g, binarize)
    make_clustersYN(self, f'eigV{which}', val, on_g, backend=backend)


def make_connected_component_by_edge(
    self: "SignedGraph",
    edge_attr='weight',
    value=1,
    on_g: str = SG_REPR,
):
    connected_components = nx.connected_components
    value_edges = [
        (u, v) for u, v, attrs in self.gr[on_g].edges(data=True) if attrs.get(edge_attr) == value
    ]
    if not value_edges:
        raise ValueError("No positive edges found in the graph.")
    G_pos = self.gr[on_g].edge_subgraph(value_edges).copy()
    all_connected_components = list(connected_components(G_pos))
    if not all_connected_components:
        raise ValueError("No connected components found with positive edges.")
    self.largest_cc = max(all_connected_components, key=len)
    self.largest_cc_subgraph = G_pos.subgraph(self.largest_cc).copy()


def handle_no_clust(self: "SignedGraph", NoClust: int) -> int | None:
    try:
        if self.numClustersBig == 0:
            logger.error("No clusters found. (SG_ERRMSG_ZEROCLUST)")
            raise NoClustError("No clusters found. (SG_ERRMSG_ZEROCLUST)")
        elif NoClust > self.numClustersBig:
            logger.error("Requested clusters exceed available. (SG_ERRMSG_NOCLUST_GTR_NCB)")
            raise NoClustError("Requested clusters exceed available. (SG_ERRMSG_NOCLUST_GTR_NCB)")
    except NoClustError as exception:
        logger.error(str(exception))
        if "GTR_NCB" in str(exception):
            msg = f"Returning the number of biggest clusters: {len(self.biggestClSet)}"
            logger.error(msg)
            return len(self.biggestClSet)
        elif "ZEROCLUST" in str(exception):
            logger.error("No clusters found, returning None.")
            return None
    return NoClust

