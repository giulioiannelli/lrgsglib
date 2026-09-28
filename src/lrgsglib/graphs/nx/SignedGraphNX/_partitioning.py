import logging
import warnings
from typing import Optional, TYPE_CHECKING
from numbers import Number

import numpy as np
import networkx as nx

from ....config.const import SG_REPR, SG_WARNMSG_NOCLUST
from ....config.errwar import NoClustError, SignedGraphWarning
from ....utils.tools.ConditionalPartitioning import *

logger = logging.getLogger(__name__)

if TYPE_CHECKING:
    from .SignedGraphNX import SignedGraphNX


def get_subgraph_from_nodes(
    self: "SignedGraphNX",
    list_of_nodes, 
    on_g: str = SG_REPR
):
    return self.gr[on_g].subgraph(list_of_nodes)


def get_nodes_subgraph_by_kv(
        self: "SignedGraphNX",
        k, 
        val: ConditionalPartitioningInput, 
        on_g: str = SG_REPR
):
    """
    Optimized version: partition nodes first, then create subgraphs.
    This avoids copying the entire graph and removing nodes one by 
    one.
    
    Parameters
    ----------
    k : str
        Node attribute key to partition on.
    val : ConditionalPartitioningInput
        Conditional partitioning specification. Can be a ConditionalPartitioning object,
        a number, a string (e.g., '>0'), or a callable.
    on_g : str, default SG_REPR
        Graph representation to use.
        
    Returns
    -------
    tuple[nx.Graph, nx.Graph]
        Two subgraphs: (G_yes, G_no) where G_yes contains nodes where the condition
        is False, and G_no contains nodes where the condition is True.
    """
    # Normalize val to ConditionalPartitioning
    val = _normalize_conditional_partitioning(val)
    
    G = self.gr[on_g]
    predicate = val.cond_func
    
    # Ensure predicate is callable
    if not callable(predicate):
        # If predicate is not callable, create a lambda that compares against it
        predicate_value = predicate
        predicate = lambda x: x == predicate_value
    
    # Partition nodes into two sets based on predicate
    nodes_yes = []  # nodes where predicate is True
    nodes_no = []  # nodes where predicate is False
    
    for node, v in G.nodes(data=k):
        if predicate(v):
            nodes_yes.append(node)
        else:
            nodes_no.append(node)
    
    # Create subgraphs using .subgraph() which is much faster
    # Note: subgraph() returns a view, use .copy() if you need 
    # independent graphs
    G_yes = G.subgraph(nodes_yes).copy()
    G_no = G.subgraph(nodes_no).copy()
    
    return G_yes, G_no


def get_eigV_cluster_sizes(
        self: "SignedGraphNX",
        which: int = 0, 
        binarize: bool = True, 
        val: ConditionalPartitioningInput = 1,
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
    val : ConditionalPartitioningInput, default 1
        Conditional partitioning specification. Can be a ConditionalPartitioning object,
        a number, a string (e.g., '>0'), or a callable. Defaults to +1.
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
    
    # Normalize val to ConditionalPartitioning
    val = _normalize_conditional_partitioning(val)
    
    if not self.clustersY:
        make_clustersYN(self, f"eigV{which}", val, on_g, backend)
    
    # Ensure clustersY exists and has valid data
    cl_len = sorted(map(len, self.clustersY), reverse=True)
    return cl_len

def get_cluster_distribution(
        self: "SignedGraphNX",
        which: int = 0, 
        on_g: str = SG_REPR,
        binarize: bool = True,
        val: ConditionalPartitioningInput = 1
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
    val : ConditionalPartitioningInput, default 1
        Conditional partitioning specification. Can be a ConditionalPartitioning object,
        a number, a string (e.g., '>0'), or a callable. Defaults to +1.
        
    Returns
    -------
    dict
        Dictionary mapping cluster sizes to their frequency.
    """
    val = _normalize_conditional_partitioning(val)
    
    cl_len = get_eigV_cluster_sizes(
        self, which, binarize, val, on_g
    )
    dictdist_cluster_sizes = {
        size: cl_len.count(size) for size in set(cl_len)
    }
    return dictdist_cluster_sizes


def get_ferroAntiferro_regions(
        self: "SignedGraphNX",
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
    self: "SignedGraphNX",
    k,
    val: ConditionalPartitioningInput,
    on_g: str = SG_REPR
):
    """
    Create Y/N graph partitions based on a conditional partitioning.
    
    Parameters
    ----------
    k : str
        Node attribute key to partition on.
    val : ConditionalPartitioningInput
        Conditional partitioning specification.
        a number, a string (e.g., '>0'), or a callable.
    on_g : str, default SG_REPR
        Graph representation to use.
    """
    # Normalize val to ConditionalPartitioning
    val = _normalize_conditional_partitioning(val)
    
    subgraph_result = self.get_nodes_subgraph_by_kv(k, val, on_g)
    self.graph_clustering_utility[k][val.key][on_g] = subgraph_result


def make_clustersYN(
    self: "SignedGraphNX",
    k,
    val: ConditionalPartitioningInput,
    on_g: str = SG_REPR,
    backend: str = 'nx',
):
    """
    Populate Y/N clusters, optionally using SciPy or CuPy connected components.
    
    Parameters
    ----------
    k : str
        Node attribute key to partition on.
    val : ConditionalPartitioningInput
        Conditional partitioning specification.
        a number, a string (e.g., '>0'), or a callable.
    on_g : str, default SG_REPR
        Graph representation to use.
    backend : str, default 'nx'
        Backend for connected components: 'nx', 'scipy', or 'cupy'.
    """
    # Normalize val to ConditionalPartitioning
    val = _normalize_conditional_partitioning(val)
    
    # The Y/N partition is cached per (k, val, on_g). For a fresh graph the
    # cache is empty on first touch (NestedDict autovivifies to {}), which is
    # normal control flow -- not an error -- so build it quietly. Only a
    # genuinely malformed entry (not a 2-tuple) falls through to the unpack
    # below, where it raises a real, visible error instead of being swallowed.
    cached = self.graph_clustering_utility[k][val.key][on_g]
    if not (isinstance(cached, (tuple, list)) and len(cached) == 2):
        make_graphYN(self, k, val, on_g)
        cached = self.graph_clustering_utility[k][val.key][on_g]
    graphY, graphN = cached

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
            from cupyx.scipy.sparse import csr_matrix as cx_csr
            from scipy.sparse.csgraph import connected_components as cs_cc

            # CuPy doesn't have connected_components, so convert to CPU
            n_comp, labels = cs_cc(adj, directed=False, return_labels=True)
            labels = np.asarray(labels, dtype=np.int64)
        else:
            raise ValueError("backend must be one of {'nx', 'scipy', 'cupy'}")

        clusters = [set() for _ in range(int(n_comp))]
        for node, lab in zip(nodes, labels):
            clusters[int(lab)].add(node)
        return clusters

    clustersY: list[set] = sorted(_components(graphY), key=len, reverse=True)
    clustersN: list[set] = _components(graphN)

    self.clustersY = clustersY
    self.clustersN = clustersN

    self.numClustersY = len(clustersY)
    self.numClustersN = len(clustersN)

    if not clustersY:
        warnings.warn(SG_WARNMSG_NOCLUST, SignedGraphWarning)
        self.biggestClSet = []
        self.numClustersBig = 0
        self.gc = None
        return

    self.biggestClSet = clustersY
    self.numClustersBig = len(clustersY)
    self.gc = clustersY[0]

    # Store graphs in the clustering utility structure
    if k not in self.graph_clustering_utility:
        self.graph_clustering_utility[k] = {}
    if val.key not in self.graph_clustering_utility[k]:
        self.graph_clustering_utility[k][val.key] = {}
    self.graph_clustering_utility[k][val.key][on_g] = (graphY, graphN)


def make_eigVclustersYN(
    self: "SignedGraphNX",
    val: ConditionalPartitioningInput,
    which: int = 0,
    on_g: str = SG_REPR,
    binarize: bool = True,
    backend: str = 'nx',
):
    # Normalize val to ConditionalPartitioning
    val = _normalize_conditional_partitioning(val)
    self.load_eigV_on_graph(which, on_g, binarize)
    self.make_clustersYN(f'eigV{which}', val, on_g, backend=backend)


def make_connected_component_by_edge(
    self: "SignedGraphNX",
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


def handle_no_clust(self: "SignedGraphNX", NoClust: int) -> int | None:
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



# ---------------------------------------------------------------------------
# k-way partitions: signed spectral clustering (Kunegis et al. 2010) and the
# heat-kernel / synchronization partitions at a diffusion scale τ.
# ---------------------------------------------------------------------------
def partition_from_features(
        self: "SignedGraphNX",
        F: np.ndarray,
        n_clusters: int,
        method: str = "kmeans",
        seed: Optional[int] = None,
        linkage_method: str = "ward",
) -> np.ndarray:
    """
    Label the nodes from a node-feature matrix ``F`` (one row per node):
    ``'kmeans'`` on the rows, or ``'linkage'``: Ward linkage on the rows
    (Euclidean, same geometry as k-means) when ``linkage_method='ward'``,
    otherwise ``linkage_method`` linkage on the distance ``1 − F`` (for a
    square similarity ``F`` in [−1, 1]). Used by ``get_spectral_partition``
    / ``get_kernel_partition``; call it directly for any other embedding
    (e.g. a SPONGE or adjacency basis).

    Returns
    -------
    np.ndarray, shape (N,) of int labels
    """
    if method == "kmeans":
        import warnings
        from sklearn.cluster import KMeans
        from sklearn.exceptions import ConvergenceWarning
        km = KMeans(n_clusters=n_clusters, n_init=10, random_state=seed)
        with warnings.catch_warnings():
            # collapsed features (e.g. a kernel at τ → ∞) give fewer distinct
            # points than clusters; the degenerate labels are the answer.
            warnings.simplefilter("ignore", ConvergenceWarning)
            return km.fit_predict(F)
    if method == "linkage":
        from scipy.cluster.hierarchy import fcluster, linkage
        from scipy.spatial.distance import squareform
        if linkage_method == "ward":
            Z = linkage(F, method="ward")
        else:
            D = 1.0 - F
            D = 0.5 * (D + D.T)
            np.fill_diagonal(D, 0.0)
            Z = linkage(squareform(D, checks=False), method=linkage_method)
        return fcluster(Z, n_clusters, criterion="maxclust") - 1
    raise ValueError(f"method must be 'kmeans' or 'linkage', got {method!r}")


def get_spectral_partition(
        self: "SignedGraphNX",
        n_clusters: int,
        n_modes: Optional[int] = None,
        laplacian_type: str = "sym",
        seed: Optional[int] = None,
) -> np.ndarray:
    """
    Signed spectral clustering: ``n_modes`` smallest eigenvectors of the signed
    Laplacian (``'signed'``: L̄ = D̄ − A; ``'sym'``: L̄_sym, Kunegis et al. SDM
    2010 §3) as node features, then k-means into ``n_clusters``.

    Computed locally: the object's cached full spectrum is left untouched.

    Returns
    -------
    np.ndarray, shape (N,) of int labels
    """
    from scipy.linalg import eigh
    M, _ = self._laplacian_operator(laplacian_type)
    M = M.toarray() if hasattr(M, "toarray") else np.asarray(M)
    m = n_clusters if n_modes is None else int(n_modes)
    _, V = eigh(M, subset_by_index=[0, m - 1])
    return self.partition_from_features(V, n_clusters, "kmeans", seed)


def get_kernel_partition(
        self: "SignedGraphNX",
        tau: float,
        n_clusters: int,
        method: str = "kmeans",
        mod_sign: bool = False,
        n_modes: Optional[int] = None,
        seed: Optional[int] = None,
        linkage_method: str = "ward",
) -> np.ndarray:
    """
    Partition at diffusion scale τ from the synchronization correlation
    ``ρ = K_τ / sqrt(diag ⊗ diag)`` (see ``get_kernel_correlation``).

    Two notions of community:
      * ``mod_sign=False`` — phase groups: nodes with the same signed profile
        ``ρ_i·`` (in phase with the same nodes, anti-phase with the same nodes);
      * ``mod_sign=True``  — synchronization groups: same ``|ρ_i·|`` profile,
        the sign relation between groups being read off afterwards.

    ``method='kmeans'`` clusters the rows of ρ (or |ρ|); ``method='linkage'``
    builds ``linkage_method`` linkage on the distance ``1 − ρ`` (or ``1 − |ρ|``)
    and cuts it into ``n_clusters``.

    Returns
    -------
    np.ndarray, shape (N,) of int labels
    """
    F = self.get_kernel_correlation(tau, n_modes=n_modes)
    if mod_sign:
        F = np.abs(F)
    return self.partition_from_features(F, n_clusters, method, seed, linkage_method)


def get_coarse_signed_adjacency(
        self: "SignedGraphNX",
        labels: np.ndarray,
        on_g: str = SG_REPR,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Positive and negative edge counts between the groups of ``labels``.

    Returns ``(W_plus, W_minus)``, each ``(n_groups, n_groups)``; the diagonal
    counts edges inside a group. Their difference is the coarse signed graph,
    the sign relation between communities.
    """
    labels = np.asarray(labels)
    groups = np.unique(labels)
    idx = {g: i for i, g in enumerate(groups)}
    Wp = np.zeros((len(groups), len(groups)))
    Wm = np.zeros_like(Wp)
    nodes = list(self.gr[on_g].nodes())
    pos = {u: k for k, u in enumerate(nodes)}
    for u, v, w in self.gr[on_g].edges(data="weight", default=1.0):
        a, b = idx[labels[pos[u]]], idx[labels[pos[v]]]
        W = Wp if w > 0 else Wm
        W[a, b] += 1
        if a != b:
            W[b, a] += 1
    return Wp, Wm


def score_partition(
        self: "SignedGraphNX",
        labels: np.ndarray,
        truth: Optional[np.ndarray] = None,
) -> dict:
    """
    Agreement of ``labels`` with ``truth`` (planted labels by default, from
    ``get_planted_labels`` when the graph has them): ``{'nmi', 'ari'}``.
    """
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    if truth is None:
        if not hasattr(self, "get_planted_labels"):
            raise ValueError("no planted labels on this graph: pass truth=")
        truth = self.get_planted_labels()
    return dict(nmi=normalized_mutual_info_score(truth, labels),
                ari=adjusted_rand_score(truth, labels))
