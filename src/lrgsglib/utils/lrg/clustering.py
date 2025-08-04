from ...shared import *
from ...nx_patches.funcs import *
from .spectral import *
from .infocomm import *

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
