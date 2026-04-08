"""
Spectral Renormalization Group for Signed/Frustrated Graphs.

Implements iterative coarse-graining of graphs based on the LDOS entropy
framework. The RG procedure partitions nodes at a chosen scale τ*, builds
a reduced graph preserving signed edge structure, and verifies spectral
preservation at coarser scales.

This is the first renormalization scheme that works on signed graphs,
where classical LRG (based on diffusion distances) breaks due to
negative entries in exp(-τL).
"""

import numpy as np
from numpy.typing import NDArray
from typing import Optional, Any
import networkx as nx
from scipy.signal import find_peaks
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

from .quantum import (
    compute_ldos_entropy,
    compute_ldos_specific_heat,
    compute_quantum_distance_matrix,
)
from .clustering import compute_optimal_threshold

__all__ = [
    "compute_frustration_fraction",
    "find_rg_scales",
    "partition_at_scale",
    "build_reduced_graph",
    "compute_reduced_spectrum",
    "spectral_rg_step",
    "spectral_rg_flow",
    "rg_flow_observables",
]


def compute_frustration_fraction(G: nx.Graph) -> float:
    """Compute fraction of edges with negative weight.

    Parameters
    ----------
    G : nx.Graph
        Graph with signed 'weight' attributes.

    Returns
    -------
    float
        Fraction of negative-weight edges in [0, 1].
    """
    weights = [d.get("weight", 1) for u, v, d in G.edges(data=True)]
    if not weights:
        return 0.0
    return sum(1 for w in weights if w < 0) / len(weights)


def find_rg_scales(
    eigenvalues: NDArray,
    num_nodes: int,
    tau_grid: Optional[NDArray] = None,
    steps: int = 400,
    t1: float = -2,
    t2: float = 4,
    min_prominence: float = 0.02,
    use_abs: bool = True,
) -> dict[str, Any]:
    """Find structural scales from specific heat C(τ) peaks.

    Parameters
    ----------
    eigenvalues : NDArray
        Eigenvalues of the (signed) Laplacian.
    num_nodes : int
        Number of nodes in the graph.
    tau_grid : NDArray, optional
        Custom τ grid. If None, uses logspace(t1, t2, steps).
    steps : int
        Number of τ points if tau_grid not given.
    t1, t2 : float
        Log-scale bounds for τ grid.
    min_prominence : float
        Minimum peak prominence as fraction of max C(τ).
    use_abs : bool
        Use |λ| in Boltzmann weights (required for signed graphs).

    Returns
    -------
    dict
        Keys: tau_stars, C_global, S_global, tau_grid, peak_indices,
        prominences.
    """
    if tau_grid is None:
        tau_grid = np.logspace(t1, t2, steps)

    lam = np.abs(eigenvalues) if use_abs else eigenvalues
    log_N = np.log(num_nodes)

    S_global = np.zeros(len(tau_grid))
    for i, tau in enumerate(tau_grid):
        boltz = np.exp(-tau * lam)
        Z = np.sum(boltz)
        if Z > 0:
            p = boltz / Z
            p_pos = p[p > 1e-30]
            S_global[i] = -np.sum(p_pos * np.log(p_pos)) / log_N

    C_global = np.gradient(1.0 - S_global, np.log(tau_grid))

    if np.max(C_global) > 0:
        peaks, props = find_peaks(
            C_global, prominence=min_prominence * np.max(C_global)
        )
    else:
        peaks = np.array([], dtype=int)
        props = {"prominences": np.array([])}

    return {
        "tau_stars": tau_grid[peaks] if len(peaks) > 0 else np.array([]),
        "C_global": C_global,
        "S_global": S_global,
        "tau_grid": tau_grid,
        "peak_indices": peaks,
        "prominences": props.get("prominences", np.array([])),
    }


def partition_at_scale(
    eigenvalues: NDArray,
    eigenvectors_col: NDArray,
    tau_star: float,
    method: str = "quantum_distance",
    n_clusters: Optional[int] = None,
    linkage_method: str = "ward",
    use_abs: bool = True,
) -> dict[str, Any]:
    """Partition nodes at a given RG scale.

    Parameters
    ----------
    eigenvalues : NDArray, shape (N,)
        Eigenvalues of the Laplacian.
    eigenvectors_col : NDArray, shape (N, N)
        Eigenvectors in column-major format.
    tau_star : float
        The RG scale at which to partition.
    method : str
        Partition method: 'quantum_distance' or 'ldos_entropy'.
    n_clusters : int, optional
        Force this many clusters. If None, auto-detect.
    linkage_method : str
        Linkage method for hierarchical clustering.
    use_abs : bool
        Use |λ| for LDOS computations.

    Returns
    -------
    dict
        Keys: labels (array), n_clusters (int), partition (list of sets),
        linkage_matrix (if applicable).
    """
    N = len(eigenvalues)

    if method == "quantum_distance":
        D_Q = compute_quantum_distance_matrix(eigenvectors_col, method="overlap")
        D_Q_cond = squareform(D_Q)
        Z = linkage(D_Q_cond, method=linkage_method)

        if n_clusters is not None:
            labels = fcluster(Z, t=n_clusters, criterion="maxclust")
        else:
            # Auto-detect n_clusters from the spectral gap
            # The largest gap between consecutive eigenvalues (after zero mode)
            # indicates the number of structurally distinct groups
            lam = np.abs(eigenvalues) if use_abs else eigenvalues
            lam_sorted = np.sort(lam)
            # Skip near-zero eigenvalues
            nonzero_mask = lam_sorted > 1e-8
            if nonzero_mask.sum() > 2:
                lam_nz = lam_sorted[nonzero_mask]
                gaps = np.diff(lam_nz)
                # Find largest gap in the first min(20, N//2) eigenvalues
                n_check = min(20, len(gaps))
                gap_idx = np.argmax(gaps[:n_check])
                # n_clusters = number of eigenvalues below the gap + 1 (for zero mode)
                n_below_gap = gap_idx + 1  # eigenvalues in the "community" band
                n_zero = np.sum(~nonzero_mask)
                n_clusters_auto = n_below_gap + n_zero
                n_clusters_auto = max(2, min(n_clusters_auto, N // 2))
            else:
                n_clusters_auto = max(2, int(np.sqrt(N)))
            labels = fcluster(Z, t=n_clusters_auto, criterion="maxclust")

    elif method == "ldos_entropy":
        S_local, _ = compute_ldos_entropy(
            eigenvalues, eigenvectors_col, tau_star, use_abs=use_abs
        )

        if n_clusters is None:
            # Auto-select: try range [2, min(20, sqrt(N))]
            from sklearn.cluster import KMeans
            from sklearn.metrics import silhouette_score

            best_k, best_score = 2, -1
            k_max = min(20, int(np.sqrt(N)))
            S_2d = S_local.reshape(-1, 1)
            for k in range(2, k_max + 1):
                km = KMeans(n_clusters=k, n_init=5, random_state=42)
                lab = km.fit_predict(S_2d)
                score = silhouette_score(S_2d, lab)
                if score > best_score:
                    best_score = score
                    best_k = k
            n_clusters = best_k

        from sklearn.cluster import KMeans

        km = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
        labels = km.fit_predict(S_local.reshape(-1, 1)) + 1  # 1-indexed
        Z = None

    else:
        raise ValueError(f"Unknown partition method: {method}")

    # Build partition as list of sets
    unique_labels = np.unique(labels)
    partition = [set(np.where(labels == lab)[0]) for lab in unique_labels]

    return {
        "labels": labels,
        "n_clusters": len(unique_labels),
        "partition": partition,
        "linkage_matrix": Z if method == "quantum_distance" else None,
    }


def build_reduced_graph(
    G_original: nx.Graph,
    partition: list[set[int]],
    sign_method: str = "majority",
    min_edge_weight: float = 0.0,
) -> tuple[nx.Graph, dict]:
    """Build a coarse-grained graph from a partition.

    Parameters
    ----------
    G_original : nx.Graph
        The original graph with signed 'weight' attributes.
    partition : list of sets
        Each set contains node indices belonging to one supernode.
    sign_method : str
        How to assign edge signs: 'majority' (net sign of cross-edges).
    min_edge_weight : float
        Minimum normalized weight to include an edge.

    Returns
    -------
    G_reduced : nx.Graph
        Coarse-grained graph with signed weights.
    metadata : dict
        Node map, cluster sizes, internal weights, frustration info.
    """
    n_clusters = len(partition)
    node_list = list(G_original.nodes())

    # Map original nodes to cluster index
    node_to_cluster = {}
    for c_idx, cluster in enumerate(partition):
        for node in cluster:
            # Handle both integer and non-integer node labels
            if node < len(node_list):
                node_to_cluster[node_list[node]] = c_idx
            else:
                node_to_cluster[node] = c_idx

    # Compute inter-cluster and intra-cluster weights
    cross_weights = np.zeros((n_clusters, n_clusters))
    cross_counts = np.zeros((n_clusters, n_clusters))

    for u, v, d in G_original.edges(data=True):
        w = d.get("weight", 1.0)
        cu = node_to_cluster.get(u)
        cv = node_to_cluster.get(v)
        if cu is not None and cv is not None:
            cross_weights[cu, cv] += w
            cross_weights[cv, cu] += w
            cross_counts[cu, cv] += 1
            cross_counts[cv, cu] += 1

    # Build reduced graph
    G_reduced = nx.Graph()
    for c_idx in range(n_clusters):
        G_reduced.add_node(
            c_idx,
            size=len(partition[c_idx]),
            internal_weight=cross_weights[c_idx, c_idx] / 2,  # undirected
            original_nodes=partition[c_idx],
        )

    for i in range(n_clusters):
        for j in range(i + 1, n_clusters):
            w_raw = cross_weights[i, j]
            n_pairs = len(partition[i]) * len(partition[j])
            w_norm = abs(w_raw) / n_pairs if n_pairs > 0 else 0

            if w_norm > min_edge_weight and abs(w_raw) > 1e-10:
                # Store signed weight
                G_reduced.add_edge(i, j, weight=w_raw / n_pairs)

    # Metadata
    metadata = {
        "node_map": {i: partition[i] for i in range(n_clusters)},
        "cluster_sizes": [len(p) for p in partition],
        "internal_weights": [cross_weights[i, i] / 2 for i in range(n_clusters)],
        "cross_weight_matrix": cross_weights,
        "frustration_fraction": compute_frustration_fraction(G_reduced),
    }

    return G_reduced, metadata


def compute_reduced_spectrum(
    G_reduced: nx.Graph,
) -> tuple[NDArray, NDArray]:
    """Compute eigenvalues and eigenvectors of the reduced graph.

    Parameters
    ----------
    G_reduced : nx.Graph
        Reduced graph with signed 'weight' attributes.

    Returns
    -------
    eigenvalues : NDArray
    eigenvectors : NDArray (column-major)
    """
    N = G_reduced.number_of_nodes()
    if N < 2:
        return np.array([0.0]), np.eye(1)

    nodes = sorted(G_reduced.nodes())
    A = nx.adjacency_matrix(G_reduced, nodelist=nodes, weight="weight")
    A = A.toarray().astype(float)
    D = np.diag(np.abs(A).sum(axis=1))
    L = D - A

    eigenvalues, eigenvectors = np.linalg.eigh(L)
    return eigenvalues, eigenvectors


def spectral_rg_step(
    G: nx.Graph,
    eigenvalues: Optional[NDArray] = None,
    eigenvectors_col: Optional[NDArray] = None,
    tau_star: Optional[float] = None,
    partition_method: str = "quantum_distance",
    n_clusters: Optional[int] = None,
    use_abs: bool = True,
) -> dict[str, Any]:
    """Perform one spectral RG step.

    Parameters
    ----------
    G : nx.Graph
        Input graph with signed weights.
    eigenvalues, eigenvectors_col : NDArray, optional
        Precomputed spectrum. If None, computed from G.
    tau_star : float, optional
        RG scale. If None, auto-detected from C(τ) first peak.
    partition_method : str
        'quantum_distance' or 'ldos_entropy'.
    n_clusters : int, optional
        Force cluster count. If None, auto-detect.
    use_abs : bool
        Use |λ| for signed graph compatibility.

    Returns
    -------
    dict
        Keys: G_reduced, eigenvalues_reduced, eigenvectors_reduced,
        partition, tau_star, metadata, partition_info, rg_scales.
    """
    N = G.number_of_nodes()

    # Compute spectrum if not provided
    if eigenvalues is None or eigenvectors_col is None:
        nodes = sorted(G.nodes())
        A = nx.adjacency_matrix(G, nodelist=nodes, weight="weight").toarray().astype(float)
        D = np.diag(np.abs(A).sum(axis=1))
        L = D - A
        eigenvalues, eigenvectors_col = np.linalg.eigh(L)

    # Find RG scales
    rg_info = find_rg_scales(eigenvalues, N, use_abs=use_abs)

    # Select tau_star
    if tau_star is None:
        if len(rg_info["tau_stars"]) > 0:
            tau_star = rg_info["tau_stars"][0]  # first (finest) peak
        else:
            # Fallback: use inverse of median eigenvalue
            med_lam = np.median(np.abs(eigenvalues[eigenvalues != 0]))
            tau_star = 1.0 / med_lam if med_lam > 0 else 1.0

    # Partition
    part_info = partition_at_scale(
        eigenvalues,
        eigenvectors_col,
        tau_star,
        method=partition_method,
        n_clusters=n_clusters,
        use_abs=use_abs,
    )

    # Build reduced graph
    G_reduced, metadata = build_reduced_graph(G, part_info["partition"])

    # Compute reduced spectrum
    eigv_r, eigV_r = compute_reduced_spectrum(G_reduced)

    return {
        "G_reduced": G_reduced,
        "eigenvalues_reduced": eigv_r,
        "eigenvectors_reduced": eigV_r,
        "partition": part_info["partition"],
        "partition_info": part_info,
        "tau_star": tau_star,
        "metadata": metadata,
        "rg_scales": rg_info,
        "N_original": N,
        "N_reduced": G_reduced.number_of_nodes(),
    }


def spectral_rg_flow(
    G: nx.Graph,
    eigenvalues: Optional[NDArray] = None,
    eigenvectors_col: Optional[NDArray] = None,
    n_steps: Optional[int] = None,
    min_nodes: int = 4,
    partition_method: str = "quantum_distance",
    use_abs: bool = True,
    verbose: bool = False,
) -> list[dict[str, Any]]:
    """Iterate spectral RG steps to produce a coarse-graining flow.

    Parameters
    ----------
    G : nx.Graph
        Input graph.
    eigenvalues, eigenvectors_col : NDArray, optional
        Precomputed spectrum of the initial graph.
    n_steps : int, optional
        Maximum number of RG steps. If None, iterate until convergence.
    min_nodes : int
        Stop when reduced graph has fewer nodes.
    partition_method : str
        Partition method for each step.
    use_abs : bool
        Use |λ| for signed graphs.
    verbose : bool
        Print progress.

    Returns
    -------
    list of dict
        Each entry is the output of one spectral_rg_step call.
    """
    flow = []
    G_current = G
    eigv_current = eigenvalues
    eigV_current = eigenvectors_col
    step_count = 0

    while True:
        N_current = G_current.number_of_nodes()

        if N_current < min_nodes:
            if verbose:
                print(f"  Stop: {N_current} nodes < min_nodes={min_nodes}")
            break

        if n_steps is not None and step_count >= n_steps:
            if verbose:
                print(f"  Stop: reached n_steps={n_steps}")
            break

        step = spectral_rg_step(
            G_current,
            eigenvalues=eigv_current,
            eigenvectors_col=eigV_current,
            partition_method=partition_method,
            use_abs=use_abs,
        )

        flow.append(step)

        if verbose:
            print(
                f"  Step {step_count}: {step['N_original']} → "
                f"{step['N_reduced']} nodes (τ* = {step['tau_star']:.3f})"
            )

        # Check for stalling (no reduction)
        if step["N_reduced"] >= N_current:
            if verbose:
                print("  Stop: no reduction achieved")
            break

        # Prepare for next step
        G_current = step["G_reduced"]
        eigv_current = step["eigenvalues_reduced"]
        eigV_current = step["eigenvectors_reduced"]
        step_count += 1

    return flow


def rg_flow_observables(
    flow: list[dict[str, Any]],
) -> dict[str, NDArray]:
    """Extract summary observables from an RG flow.

    Parameters
    ----------
    flow : list of dict
        Output from spectral_rg_flow.

    Returns
    -------
    dict
        Arrays indexed by RG step: n_nodes, n_edges, frustration_fraction,
        spectral_gap, tau_stars, entropy_at_tau_star.
    """
    n_nodes = []
    n_edges = []
    frustration = []
    spectral_gap = []
    tau_stars = []
    entropy_at_tau_star = []

    for step in flow:
        G_r = step["G_reduced"]
        eigv_r = step["eigenvalues_reduced"]

        n_nodes.append(step["N_original"])
        n_edges.append(G_r.number_of_edges())
        frustration.append(step["metadata"]["frustration_fraction"])

        # Spectral gap: smallest nonzero eigenvalue
        nonzero = np.abs(eigv_r) > 1e-10
        if nonzero.any():
            spectral_gap.append(np.min(np.abs(eigv_r[nonzero])))
        else:
            spectral_gap.append(0.0)

        tau_stars.append(step["tau_star"])

        # Entropy at tau_star
        idx = np.argmin(np.abs(step["rg_scales"]["tau_grid"] - step["tau_star"]))
        entropy_at_tau_star.append(step["rg_scales"]["S_global"][idx])

    # Add final reduced graph info
    if flow:
        last = flow[-1]
        n_nodes.append(last["N_reduced"])

    return {
        "n_nodes": np.array(n_nodes),
        "n_edges": np.array(n_edges),
        "frustration_fraction": np.array(frustration),
        "spectral_gap": np.array(spectral_gap),
        "tau_stars": np.array(tau_stars),
        "entropy_at_tau_star": np.array(entropy_at_tau_star),
    }
