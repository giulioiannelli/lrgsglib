from typing import Any, Dict, List, Optional, Tuple

import networkx as nx
import numpy as np
from numpy.typing import NDArray
from scipy.cluster.hierarchy import dendrogram, linkage

from ...graphs.nx.funcs.spectral import signed_laplacian_matrix
from .infocomm import lapl_dists

__all__ = [
    "MakeLinkageMatrix",
    "compute_normalized_linkage",
    "compute_optimal_threshold",
    "circular_layout_by_cluster",
    "log_dendrogram",
    "dendrogram_leaf_node_colors",
]

def MakeLinkageMatrix(G: nx.Graph, tau: float = 1e-2, is_signed: bool = False, method: str = "ward") -> Tuple[np.ndarray, np.ndarray]:
    """
    Constructs a hierarchical clustering linkage matrix for a graph.

    Parameters
    ----------
    G : networkx.Graph
        The input graph.
    tau : float, optional
        A threshold parameter for distance computation (default is 1e-2).
    is_signed : bool, optional
        Whether the graph is signed (default is False).
    method : str, optional
        Linkage method to use (default is "ward").

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        The hierarchical clustering linkage matrix and the weights associated with the graph's Laplacian spectrum.
    """
    L = signed_laplacian_matrix(G)
    dists = lapl_dists(L, tau, is_signed)
    linkage_matrix1 = linkage(dists, method=method)
    tmax = linkage_matrix1[::, 2][-1]
    linkage_matrix = linkage(dists / tmax, method=method)
    return linkage_matrix


def compute_normalized_linkage(dists: np.ndarray, G: nx.Graph, method: str = "average", labelList: str = "names") -> Tuple[np.ndarray, List, float]:
    """
    Computes a normalized hierarchical clustering linkage matrix.

    Parameters
    ----------
    dists : np.ndarray
        Condensed distance matrix (e.g., output from scipy.spatial.distance.pdist).
    G : networkx.Graph
        Graph whose nodes correspond to the distances in `dists`.
    method : str, optional
        Linkage method to use (default is "average").
    labelList : str, optional
        Specifies the type of labels to generate for nodes. Options are "names" (default) or "numbers".

    Returns
    -------
    Tuple[np.ndarray, List, float]
        The normalized hierarchical clustering linkage matrix, the list of node labels, and the normalization factor used.
    """
    # Initial linkage computation.
    linkage_matrix = linkage(dists, method)
    
    # Create labels for nodes.
    if labelList == "names":
        label_list = list(G.nodes())
    elif labelList == "numbers":
        label_list = [i + 1 for i in range(len(G.nodes()))]
    
    # Determine tmax as the maximum merge distance plus 1% of that distance.
    max_distance = linkage_matrix[-1, 2]
    tmax = max_distance + 0.01 * max_distance
    
    # Recompute the linkage matrix using normalized distances.
    linkage_matrix = linkage(dists / tmax, method)
    
    return linkage_matrix, label_list, tmax


def compute_optimal_threshold(linkage_matrix: np.ndarray, scaling_factor: float = 0.9) -> Tuple[float, float, np.ndarray, int]:
    """
    Compute the optimal flat clustering threshold from a linkage matrix
    using the partition stability index.

    Parameters
    ----------
    linkage_matrix : np.ndarray
        The linkage matrix from hierarchical clustering. It is assumed that 
        the third column contains the merge distances.
    scaling_factor : float, optional
        A factor to scale the optimal threshold (default is 0.9).

    Returns
    -------
    Tuple[float, float, np.ndarray, int]
        The computed threshold for flat clustering, the optimal threshold derived from the dendrogram gap analysis,
        the array of computed stability indices for each branch, and the index corresponding to the branch with the maximum stability index.
    """
    # Extract merge distances from the linkage matrix
    dendro_thresholds = linkage_matrix[:, 2]
    # Reverse order so that the first element corresponds to the initial split
    D_values = dendro_thresholds[::-1]
    
    # Compute the normalization constant N using the first and last thresholds
    N = 1 / (np.log10(D_values[0]) - np.log10(D_values[-1]))
    
    # Compute the stability index for each dendrogram gap
    stability_indices = []
    for i in range(len(D_values) - 1):
        sigma = N * (np.log10(D_values[i]) - np.log10(D_values[i+1]))
        stability_indices.append(sigma)
    stability_indices = np.array(stability_indices)
    
    # Identify the branch with the highest stability index
    optimal_branch_index = np.argmax(stability_indices)
    # The optimal threshold is D_(n+1) corresponding to that branch
    optimal_threshold = D_values[optimal_branch_index + 1]
    
    # Apply a scaling factor to determine the final flat clustering threshold
    FlatClusteringTh = optimal_threshold * scaling_factor

    return FlatClusteringTh, optimal_threshold, stability_indices, optimal_branch_index




def log_dendrogram(
    Z: NDArray,
    ax: Any = None,
    color_threshold: Optional[float] = None,
    above_threshold_color: str = "#888888",
    ylim_rel_pad: float = 0.0,
    normalize: bool = True,
    draw_cut: bool = True,
    cut_color: str = "r",
    cut_linestyle: str = "--",
    set_ylabel: bool = True,
    no_labels: bool = True,
    **kwargs: Any,
) -> Tuple[dict, float, float, float]:
    """Plot a hierarchical-clustering dendrogram with a log-scaled y-axis.

    ``scipy.cluster.hierarchy.dendrogram`` draws leaf bars starting at
    ``y = 0``, which becomes ``-inf`` on a log scale — producing very long
    leaf lines and forcing the automatic y-limit to cover an empty range.
    This helper:

    1. Normalizes the linkage heights by ``D_max = max(Z[:, 2])`` so
       the y-axis runs in ``(0, 1]`` (disable with ``normalize=False``).
    2. Post-processes the line collections drawn by
       ``scipy.cluster.hierarchy.dendrogram`` so that every segment whose
       y-coordinate falls below the smallest positive merge height is
       clamped to that value — eliminating the "long leaf line" artefact.
    3. Switches the axis to log scale and sets a tight ``ylim`` derived
       from the linkage matrix itself, padded by ``ylim_log_pad`` decades
       on each side.
    4. Optionally draws a horizontal cut line at ``color_threshold``
       (translated to the normalized scale).

    Parameters
    ----------
    Z : NDArray
        Linkage matrix as returned by ``scipy.cluster.hierarchy.linkage``.
    ax : matplotlib.axes.Axes, optional
        Axis on which to draw. If ``None`` the current axis is used.
    color_threshold : float, optional
        Merge-height threshold (in the **raw** un-normalized scale of
        ``Z``) at which to colour branches for cluster separation.
    above_threshold_color : str, optional
        Colour used for branches above the threshold (default grey
        ``#888888``).
    ylim_rel_pad : float, optional
        Padding added above and below the linkage-height range,
        expressed as a **fraction of the tree's log span**
        (default ``0.08`` = 8 %). This keeps the empty margin tight
        regardless of whether the merge heights span a fraction of a
        decade (as in lattices) or several decades (as in hierarchical
        networks).
    normalize : bool, optional
        If ``True`` (default) divide all merge heights and the
        ``color_threshold`` by ``D_max = Z[:, 2].max()`` so the y-axis
        runs in ``(0, 1]``.
    draw_cut : bool, optional
        If ``True`` (default) draw a horizontal line at the
        ``color_threshold`` (normalized if applicable).
    cut_color, cut_linestyle : str, optional
        Style of the cut line.
    set_ylabel : bool, optional
        If ``True`` (default) label the y-axis as
        :math:`\\mathcal D/\\mathcal D_{\\max}` (normalized case) or
        :math:`\\mathcal D` (raw case).
    no_labels : bool, optional
        Whether to suppress the default leaf labels (default ``True``).
    **kwargs : Any
        Additional keyword arguments forwarded to
        ``scipy.cluster.hierarchy.dendrogram``.

    Returns
    -------
    ddata : dict
        The dictionary returned by ``scipy.cluster.hierarchy.dendrogram``.
    t_min : float
        Smallest positive merge height used for the lower ``ylim``
        (in normalized units if ``normalize=True``).
    t_max : float
        Largest merge height used for the upper ``ylim``
        (``1.0`` if ``normalize=True``).
    D_max : float
        The original (un-normalized) maximum merge height
        ``Z[:, 2].max()``.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        ax = plt.gca()

    raw_heights = np.asarray(Z[:, 2], dtype=float)
    D_max = float(raw_heights.max())
    if normalize and D_max <= 0:
        raise ValueError("Cannot normalize a linkage with non-positive D_max.")

    Z_plot = np.array(Z, dtype=float, copy=True)
    if normalize:
        Z_plot[:, 2] = Z_plot[:, 2] / D_max
        color_threshold_plot = (
            color_threshold / D_max if color_threshold is not None else None
        )
    else:
        color_threshold_plot = color_threshold

    heights_plot = Z_plot[:, 2]
    positive_heights = heights_plot[heights_plot > 0]
    if positive_heights.size == 0:
        raise ValueError(
            "Linkage matrix has no positive merge heights; cannot draw "
            "a log-scaled dendrogram."
        )
    t_min = float(positive_heights.min())
    t_max = float(heights_plot.max())

    ddata = dendrogram(
        Z_plot,
        ax=ax,
        color_threshold=color_threshold_plot,
        above_threshold_color=above_threshold_color,
        no_labels=no_labels,
        **kwargs,
    )

    # Clamp any line-segment y-coordinate below ``t_min`` to ``t_min``.
    for collection in ax.collections:
        if not hasattr(collection, "get_segments"):
            continue
        segments = collection.get_segments()
        new_segments = []
        for seg in segments:
            seg = np.array(seg, dtype=float)
            seg[seg[:, 1] < t_min, 1] = t_min
            new_segments.append(seg)
        collection.set_segments(new_segments)

    ax.set_yscale("log")
    # Pad in log space by a fraction of the tree's own log span, so a
    # compact tree (lattice) gets proportionally small padding and a
    # wide tree (hierarchical) gets proportionally larger padding.
    log_span = np.log10(t_max) - np.log10(t_min)
    # Protect against degenerate (all-equal) trees.
    log_span = max(log_span, 1e-12)
    log_pad = ylim_rel_pad * log_span
    lo = t_min * (10.0 ** (-log_pad))
    hi = t_max * (10.0 ** (+log_pad))
    ax.set_ylim(lo, hi)

    if draw_cut and color_threshold_plot is not None:
        ax.axhline(
            color_threshold_plot,
            color=cut_color,
            linestyle=cut_linestyle,
            linewidth=1.0,
        )

    if set_ylabel:
        if normalize:
            ax.set_ylabel(r"$\mathcal{D}/\mathcal{D}_{\max}$")
        else:
            ax.set_ylabel(r"$\mathcal{D}$")

    return ddata, t_min, t_max, D_max


def dendrogram_leaf_node_colors(
    ddata: dict, n_nodes: int
) -> Tuple[List[str], int]:
    """Build per-node colours consistent with a scipy dendrogram.

    Extracts the leaf colours assigned by ``scipy.cluster.hierarchy.dendrogram``
    (via the ``color_threshold`` argument) and returns a list aligned with
    the original node indices ``0..n_nodes - 1``. The colours returned are
    exactly the ones used to paint the dendrogram leaves, so downstream
    plots of the same graph will match the dendrogram visually.

    Parameters
    ----------
    ddata : dict
        Dictionary returned by ``scipy.cluster.hierarchy.dendrogram``
        (or by :func:`log_dendrogram`).
    n_nodes : int
        Total number of nodes in the original data (typically ``G.number_of_nodes()``).

    Returns
    -------
    node_colors : list of str
        ``node_colors[i]`` is the dendrogram leaf colour of node ``i``.
    n_clusters : int
        Number of distinct leaf colours excluding the
        ``above_threshold_color`` (default grey). This is the effective
        number of clusters implied by the dendrogram's colour cut.
    """
    leaves = ddata["leaves"]
    colors = ddata["leaves_color_list"]
    mapping = {leaves[k]: colors[k] for k in range(len(leaves))}
    node_colors = [mapping[i] for i in range(n_nodes)]
    unique_cluster_colors = {
        c for c in colors if not (c.startswith("#8") or c == "grey" or c == "gray")
    }
    return node_colors, len(unique_cluster_colors)


def circular_layout_by_cluster(G: nx.Graph, cluster_assignment: Dict) -> Dict:
    """
    Compute a circular layout where nodes are ordered by their cluster assignment.

    Parameters
    ----------
    G : networkx.Graph
        The input graph.
    cluster_assignment : Dict
        Mapping from node to its cluster id.

    Returns
    -------
    Dict
        Dictionary mapping node to (x, y) position.
    """
    # Sort nodes by cluster id (and optionally by node id for tie-breaking)
    nodes_sorted = sorted(G.nodes(), key=lambda node: (cluster_assignment[node], node))
    n = len(nodes_sorted)
    pos = {}
    for i, node in enumerate(nodes_sorted):
        angle = 2 * np.pi * i / n
        pos[node] = (np.cos(angle), np.sin(angle))
    return pos
