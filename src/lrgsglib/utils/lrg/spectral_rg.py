"""
Spectral Renormalization Group for Signed/Frustrated Graphs.

Implements iterative coarse-graining of graphs based on the LDOS entropy
framework derived from the quantum propagator on graphs.

The key quantity is the signed diffusion distance:

    d_diff(i,j,τ)² = Σ_k exp(-2τ|λ_k|) (v_k(i) - v_k(j))²

The |λ_k| weighting is not an ad-hoc fix for signed graphs — it is the
natural consequence of the quantum propagator exp(-itL) being unitary
(|exp(-itλ)| = 1 for all λ). Time-averaging the quantum walk produces
|λ| in the Boltzmann factors, making the signed Laplacian unnecessary.

On unsigned graphs, |λ| = λ and d_diff reduces exactly to the classical
diffusion distance, recovering standard LRG. On signed/frustrated graphs,
the formalism extends naturally without blow-up.

The eigenvector structure of a single spectral decomposition encodes two
channels: topology (block/community) and geometry (frustration domains).
These are not separate computations — they are two views of the same
spectral data, separated by scale (τ) and mode structure.
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
    "negative_edge_fraction",
    "frustration_index",
    "spectral_frustration",
    "compute_frustration_fraction",  # deprecated alias
    "compute_signed_diffusion_distance",
    "compute_eigenmode_sign_distance",
    "agmon_geodesic_distance",
    "find_rg_scales",
    "partition_at_scale",
    "build_reduced_graph",
    "compute_reduced_spectrum",
    "spectral_rg_step",
    "spectral_rg_flow",
    "rg_flow_observables",
]


def negative_edge_fraction(G: nx.Graph) -> float:
    """Fraction of edges with negative weight.

    This is NOT the frustration of the graph — an all-negative bipartite
    graph has no frustrated cycles but 100% negative edges. For actual
    frustration, use :func:`frustration_index` or :func:`spectral_frustration`.

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


def frustration_index(
    G: nx.Graph, cycles: Optional[list] = None
) -> float:
    """Fraction of cycles that are frustrated (product of edge signs = -1).

    A cycle is **frustrated** iff it contains an odd number of negative
    edges. The frustration index is the fraction of cycles (from a
    supplied basis, or all triangles by default) that are frustrated.

    This is **zero** on balanced graphs: all-positive, all-negative
    bipartite lattices, Harary-balanced signed graphs, etc. It is
    **positive** only when the sign structure contains genuine
    topological obstructions.

    Parameters
    ----------
    G : nx.Graph
        Graph with signed 'weight' attributes.
    cycles : list of lists, optional
        A collection of cycles as lists of nodes forming each cycle.
        If None, uses all triangles in G. For lattices, pass the
        plaquettes.

    Returns
    -------
    float
        Fraction of frustrated cycles in [0, 1]. Returns 0 if no
        cycles are available.
    """
    if cycles is None:
        # Use triangles as default cycle basis
        cycles = []
        for u in G.nodes():
            nbrs = list(G.neighbors(u))
            for i, v in enumerate(nbrs):
                for w in nbrs[i + 1 :]:
                    if G.has_edge(v, w) and u < v and u < w:
                        cycles.append([u, v, w])

    if not cycles:
        return 0.0

    n_frustrated = 0
    for cycle in cycles:
        prod = 1.0
        k = len(cycle)
        valid = True
        for i in range(k):
            u, v = cycle[i], cycle[(i + 1) % k]
            if not G.has_edge(u, v):
                valid = False
                break
            prod *= G[u][v].get("weight", 1.0)
        if valid and prod < 0:
            n_frustrated += 1

    return n_frustrated / len(cycles)


def lattice_plaquettes(side: int, periodic: bool = True) -> list:
    """Return the list of unit plaquettes (4-cycles) of a square lattice.

    Parameters
    ----------
    side : int
        Lattice side length.
    periodic : bool
        If True, use periodic boundary conditions.

    Returns
    -------
    list of lists
        Each plaquette is a list of 4 integer node indices.
    """
    plaquettes = []
    r_max = side if periodic else side - 1
    c_max = side if periodic else side - 1
    for r in range(r_max):
        for c in range(c_max):
            r1 = (r + 1) % side
            c1 = (c + 1) % side
            plaquettes.append([
                r * side + c,
                r * side + c1,
                r1 * side + c1,
                r1 * side + c,
            ])
    return plaquettes


def spectral_frustration(G: nx.Graph) -> float:
    """Spectral frustration: smallest eigenvalue of Kunegis signed Laplacian.

    The Kunegis signed Laplacian L_s = D_|w| - A_signed is PSD. It has a
    zero eigenvalue iff the graph is structurally balanced (Harary).
    The smallest eigenvalue is a continuous measure of frustration:
    zero for balanced graphs, positive otherwise.

    Normalized by max(|λ|) to give a dimensionless quantity in [0, 1].

    Parameters
    ----------
    G : nx.Graph
        Graph with signed 'weight' attributes.

    Returns
    -------
    float
        Normalized smallest eigenvalue, zero for balanced graphs.
    """
    if G.number_of_nodes() < 2:
        return 0.0
    nodes = sorted(G.nodes())
    A = nx.adjacency_matrix(G, nodelist=nodes, weight="weight").toarray().astype(float)
    D = np.diag(np.abs(A).sum(axis=1))
    L = D - A
    eigv = np.linalg.eigvalsh(L)
    lam_max = max(abs(eigv).max(), 1e-12)
    return float(eigv[0]) / lam_max


# Deprecated alias — kept for backwards compatibility. This is NOT a
# proper frustration measure; it's just the negative edge fraction.
compute_frustration_fraction = negative_edge_fraction


def compute_signed_diffusion_distance(
    eigenvalues: NDArray,
    eigenvectors_col: NDArray,
    tau: float,
    use_abs: bool = True,
) -> NDArray:
    """Compute signed diffusion distance matrix at scale τ.

    The signed diffusion distance uses |λ_k| in the Boltzmann weights,
    making it valid for signed/frustrated graphs where classical diffusion
    distances (using exp(-τλ)) break due to negative eigenvalues.

    .. math::

        d_{\\text{diff}}(i,j,\\tau)^2 = \\sum_k e^{-2\\tau|\\lambda_k|}
        \\bigl(v_k(i) - v_k(j)\\bigr)^2

    Parameters
    ----------
    eigenvalues : NDArray, shape (N,)
        Eigenvalues of the (signed) Laplacian.
    eigenvectors_col : NDArray, shape (N_nodes, N_modes)
        Eigenvectors in column-major format: ``eigV_col[i, k] = v_k(i)``.
    tau : float
        Diffusion time / RG scale parameter.
    use_abs : bool
        If True, use |λ_k| in Boltzmann weights (required for signed graphs).

    Returns
    -------
    NDArray, shape (N_nodes, N_nodes)
        Symmetric distance matrix with zero diagonal.

    Notes
    -----
    At small τ, all modes contribute approximately equally (like d_Q).
    At large τ, only the lowest-eigenvalue modes survive, giving a
    coarse-grained distance sensitive to global structure.
    The optimal τ typically scales as τ ~ L² where L is the spatial scale
    of the feature being detected.
    """
    lam = np.abs(eigenvalues) if use_abs else eigenvalues
    phases = np.exp(-tau * lam)
    # Weighted eigenvectors: V_w[i, k] = v_k(i) * exp(-tau * |lambda_k|)
    V_w = eigenvectors_col * phases[np.newaxis, :]
    # d^2(i,j) = ||V_w[i,:] - V_w[j,:]||^2 via expansion
    diag = np.sum(V_w ** 2, axis=1)
    cross = V_w @ V_w.T
    D2 = diag[:, np.newaxis] + diag[np.newaxis, :] - 2.0 * cross
    D = np.sqrt(np.maximum(D2, 0.0))
    np.fill_diagonal(D, 0.0)
    return D


def compute_eigenmode_sign_distance(
    eigenvalues: NDArray,
    eigenvectors_col: NDArray,
    tau: float,
    cutoff: float = 1.0,
    use_abs: bool = True,
    threshold: float = 1e-12,
    smooth_sign_scale: Optional[float] = None,
) -> NDArray:
    """Compute a multiscale distance from joint eigenmode sign signatures.

    The distance is built from the cumulative sign structure of the modes
    active up to scale ``tau``:

    .. math::

        d_{\\mathrm{sign}}(i,j;\\tau)^2
        =
        \\frac{1}{W_\\tau}
        \\sum_k w_k(\\tau)
        \\frac{1 - s_k(i)s_k(j)}{2},

    where ``s_k(i) = sign(v_k(i))`` and the default cutoff weights are

    .. math::

        w_k(\\tau) = \\mathbf{1}[\\tau \\tilde\\lambda_k \\le \\mathrm{cutoff}],
        \\qquad
        \\tilde\\lambda_k =
        \\begin{cases}
            |\\lambda_k| & \\text{if use_abs=True} \\\\
            \\lambda_k   & \\text{if use_abs=False}
        \\end{cases}

    The primitive compares the *joint sign signature* of the modes that are
    still active at scale ``tau``. On open lattices this directly tracks the
    nodal-domain geometry of standing-wave eigenmodes. The resulting distance
    is symmetric, bounded, and genuinely multiscale, but it is not expected
    to be ultrametric by itself; an ultrametric can be induced afterwards via
    hierarchical clustering.

    Parameters
    ----------
    eigenvalues : NDArray, shape (N,)
        Eigenvalues of the Laplacian.
    eigenvectors_col : NDArray, shape (N_nodes, N_modes)
        Eigenvectors in column-major format: ``eigV_col[i, k] = v_k(i)``.
    tau : float
        Scale parameter.
    cutoff : float
        Keep modes satisfying ``tau * lambda_eff <= cutoff``.
    use_abs : bool
        Use ``|lambda_k|`` instead of ``lambda_k`` in the cutoff rule.
    threshold : float
        Numerical threshold below which eigenvalues are treated as zero modes
        and excluded from the sign signature.
    smooth_sign_scale : float, optional
        If provided, replace the hard sign with ``tanh(v / smooth_sign_scale)``
        to soften nodal boundaries near zeros.

    Returns
    -------
    NDArray, shape (N_nodes, N_nodes)
        Symmetric sign-geometry distance matrix with zero diagonal.
    """
    lam = np.abs(eigenvalues) if use_abs else eigenvalues
    active = (lam > threshold) & (tau * lam <= cutoff)

    if not np.any(active):
        nonzero = np.flatnonzero(lam > threshold)
        if len(nonzero) == 0:
            return np.zeros(
                (eigenvectors_col.shape[0], eigenvectors_col.shape[0]),
                dtype=float,
            )
        active[nonzero[0]] = True

    V_active = eigenvectors_col[:, active]
    if smooth_sign_scale is None:
        signs = np.sign(V_active)
    else:
        signs = np.tanh(V_active / smooth_sign_scale)

    corr = (signs @ signs.T) / signs.shape[1]
    D2 = 0.5 * (1.0 - corr)
    D2 = np.maximum(0.5 * (D2 + D2.T), 0.0)
    np.fill_diagonal(D2, 0.0)
    return np.sqrt(D2)


def agmon_geodesic_distance(
    sg: Any,
    eigenvectors_col: NDArray,
    n_modes: Optional[int] = None,
    mode_indices: Optional[Any] = None,
    eigenvalues: Optional[NDArray] = None,
    tau: Optional[float] = None,
    use_abs: bool = True,
) -> NDArray:
    """Gradient-energy geodesic (Agmon / propagator) distance on a signed graph.

    For each edge ``(a, b)`` of ``sg`` the function assigns a weight
    that is the gradient energy of a (possibly Boltzmann-weighted)
    spectral projector and then returns the all-pairs shortest-path
    distance on ``sg`` with those edge weights.

    The function supports two complementary mode-selection schemes:

    **1. Hard projector (default).** Pass ``n_modes`` or
    ``mode_indices``. The edge weight is

    .. math::

        w_{ab}^{2}
        \\;=\\;
        \\sum_{k\\in\\mathcal M}\\bigl(v_{k}(a)-v_{k}(b)\\bigr)^{2}
        \\;=\\;
        (e_{a}-e_{b})^{\\top}P_{\\mathcal M}(e_{a}-e_{b}),

    where :math:`P_{\\mathcal M}=\\sum_{k\\in\\mathcal M}|v_{k}\\rangle\\langle v_{k}|`
    is the orthogonal projector onto the selected subspace. The metric
    depends only on the projector, hence is basis-free within any
    degenerate eigenspace. Cutoffs should sit at eigenvalue-gap
    boundaries — cutting inside a multiplet leaves the projector
    basis-dependent. On inhomogeneous networks the hard projector can
    produce singleton-cluster artefacts: high-frequency modes whose
    nodal sets graze single nodes contribute equally to the metric and
    rip those nodes out of their natural community.

    **2. Soft Boltzmann / propagator weighting.** Pass ``tau`` (and
    ``eigenvalues``). The edge weight becomes

    .. math::

        w_{ab}^{2}(\\tau)
        \\;=\\;
        \\sum_{k}e^{-2\\tau|\\lambda_{k}|}\\bigl(v_{k}(a)-v_{k}(b)\\bigr)^{2}
        \\;=\\;
        (e_{a}-e_{b})^{\\top}\\,e^{-2\\tau|\\hat L|}\\,(e_{a}-e_{b}),

    i.e. the squared signed-diffusion (LRG) distance promoted to an
    edge metric. The hard projector :math:`P_{\\mathcal M}` is replaced
    by the smooth propagator :math:`e^{-2\\tau|\\hat L|}` so high-frequency
    modes are exponentially damped — this kills the singleton artefacts
    and recovers a clean nested partition on inhomogeneous networks.
    The :math:`|\\hat L|` form makes the construction valid on signed
    graphs as well; on unsigned graphs it reduces to the classical
    diffusion edge metric (:math:`D\\sim 1/\\hat\\rho_{ij}` of LRG, made
    into a geodesic).

    The two schemes can be combined: pass *both* a mode subset
    (``n_modes`` / ``mode_indices``) and ``tau``, in which case the
    Boltzmann factors are applied only to the retained modes.

    Engine-agnostic: accepts any object satisfying
    :class:`~lrgsglib.graphs.protocols.SignedGraphProtocol` (both
    ``SignedGraphNX`` and ``SignedGraphGT`` and all their subclasses).
    Only ``sg.N`` and ``sg.get_edges_with_weights()`` are used, so no
    backend-specific graph object is ever handled here.

    Parameters
    ----------
    sg : SignedGraphProtocol
        Signed graph object (any backend). Only the edge list is used;
        existing edge weights are ignored — the metric comes entirely
        from the eigenmodes.
    eigenvectors_col : NDArray, shape (N, N_modes)
        Eigenvectors in column-major format: ``eigV[i, k] = v_k(i)``.
        Row :math:`i` corresponds to node index :math:`i` (the same
        convention used by ``sg.get_edges_with_weights()``).
    n_modes : int, optional
        Use the first ``n_modes`` columns of ``eigenvectors_col``.
        Mutually exclusive with ``mode_indices``.
    mode_indices : sequence of int, optional
        Explicit list of mode column indices to include.
    eigenvalues : NDArray, optional
        Eigenvalues, aligned with the columns of ``eigenvectors_col``.
        Required when ``tau`` is given.
    tau : float, optional
        Diffusion / RG scale. If supplied, applies the Boltzmann
        factor :math:`e^{-2\\tau|\\lambda_k|}` to each retained mode.
    use_abs : bool, default True
        Use :math:`|\\lambda_k|` instead of :math:`\\lambda_k` in the
        Boltzmann factor. Required for signed graphs (where some
        eigenvalues are negative); harmless on unsigned graphs.

    Returns
    -------
    NDArray, shape (N, N)
        Symmetric geodesic-distance matrix with zero diagonal.

    Notes
    -----
    - Edge construction is :math:`O(E\\cdot n_\\text{modes})` and
      Dijkstra is :math:`O((N+E)\\log N)`.
    - The distance is a proper metric on ``sg``: symmetric,
      non-negative, and satisfies the triangle inequality by
      construction.
    - In the limit ``tau -> 0`` and using all modes, the Boltzmann
      form collapses to :math:`w_{ab}=\\sqrt 2` for every edge, so
      the geodesic reduces to :math:`\\sqrt 2` times the unweighted
      hop distance — pure topological community detection.
    """
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import shortest_path

    if (n_modes is not None) and (mode_indices is not None):
        raise ValueError(
            "agmon_geodesic_distance: pass at most one of "
            "`n_modes` or `mode_indices`."
        )
    if tau is not None and eigenvalues is None:
        raise ValueError(
            "agmon_geodesic_distance: `tau` requires `eigenvalues`."
        )
    if tau is None and n_modes is None and mode_indices is None:
        raise ValueError(
            "agmon_geodesic_distance: pass at least one of "
            "`n_modes`, `mode_indices`, or `tau` (with `eigenvalues`)."
        )

    if n_modes is not None:
        V = eigenvectors_col[:, : int(n_modes)]
        if eigenvalues is not None:
            lam_sel = np.asarray(eigenvalues)[: int(n_modes)]
        else:
            lam_sel = None
    elif mode_indices is not None:
        idx = list(mode_indices)
        V = eigenvectors_col[:, idx]
        if eigenvalues is not None:
            lam_sel = np.asarray(eigenvalues)[idx]
        else:
            lam_sel = None
    else:
        V = eigenvectors_col
        lam_sel = np.asarray(eigenvalues)

    if tau is not None:
        lam = np.abs(lam_sel) if use_abs else lam_sel
        weights = np.exp(-2.0 * float(tau) * lam)
        V = V * np.sqrt(np.maximum(weights, 0.0))[np.newaxis, :]

    N = int(sg.N)
    edges = sg.get_edges_with_weights()
    rows: list[int] = []
    cols: list[int] = []
    data: list[float] = []
    for u, v, _w in edges:
        iu = int(u)
        iv = int(v)
        w2 = float(((V[iu] - V[iv]) ** 2).sum())
        w = float(np.sqrt(max(w2, 0.0)))
        rows += [iu, iv]
        cols += [iv, iu]
        data += [w, w]
    W = csr_matrix((data, (rows, cols)), shape=(N, N))
    return shortest_path(W, method="D", directed=False)


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
    method: str = "joint",
    n_clusters: Optional[int] = None,
    linkage_method: str = "ward",
    use_abs: bool = True,
    G: Optional[nx.Graph] = None,
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
        Partition method:

        - ``'signed_diffusion'`` (default): Hierarchical clustering on
          the signed diffusion distance using the provided eigenvectors.
        - ``'topological'``: Same d_diff formula but recomputes the
          eigendecomposition from the **unsigned** Laplacian (|weights|).
          Detects block/community structure that survives frustration.
          Requires ``G`` parameter.
        - ``'joint'``: Computes both signed and topological d_diff
          embeddings and combines them into one product-space distance.
          Two nodes are close iff they are close in BOTH embeddings —
          same connectivity community AND same sign domain. This is
          the most principled partition method as it does not hide
          either information source. Requires ``G`` parameter.
        - ``'eigenmode_sign'``: Clusters nodes by low-mode sign pattern.
          Groups nodes into frustration-consistent domains (ferro or
          antiferro), keeping frustration at cluster boundaries rather
          than absorbing it. Best for frustration-preserving RG.
          ``tau_star`` controls how many modes to include via
          Boltzmann weighting.
        - ``'quantum_distance'``: Hierarchical clustering on d_Q (legacy).
        - ``'ldos_entropy'``: K-means on per-node S(τ,i).

    n_clusters : int, optional
        Force this many clusters. If None, auto-detect via silhouette scan.
    linkage_method : str
        Linkage method for hierarchical clustering.
    use_abs : bool
        Use |λ| in Boltzmann weights (required for signed graphs).
    G : nx.Graph, optional
        Required for ``method='topological'``. The graph from which to
        compute unsigned Laplacian eigenvectors.

    Returns
    -------
    dict
        Keys: labels (array), n_clusters (int), partition (list of sets),
        linkage_matrix (ndarray or None), distance_matrix (ndarray or None).
    """
    N = len(eigenvalues)

    if method == "signed_diffusion":
        D = compute_signed_diffusion_distance(
            eigenvalues, eigenvectors_col, tau_star, use_abs=use_abs
        )
        D_cond = squareform(D)
        Z = linkage(D_cond, method=linkage_method)

        if n_clusters is not None:
            labels = fcluster(Z, t=n_clusters, criterion="maxclust")
        else:
            labels = _auto_k_silhouette(D, Z, N)

    elif method == "topological":
        if G is None:
            raise ValueError("method='topological' requires the G parameter")
        nodes = sorted(G.nodes())
        A_u = nx.adjacency_matrix(G, nodelist=nodes, weight="weight")
        A_u = np.abs(A_u.toarray().astype(float))
        D_u = np.diag(A_u.sum(axis=1))
        L_u = D_u - A_u
        eigv_u, eigV_u = np.linalg.eigh(L_u)

        D = compute_signed_diffusion_distance(
            eigv_u, eigV_u, tau_star, use_abs=False
        )
        D_cond = squareform(D)
        Z = linkage(D_cond, method=linkage_method)

        if n_clusters is not None:
            labels = fcluster(Z, t=n_clusters, criterion="maxclust")
        else:
            labels = _auto_k_silhouette(D, Z, N)

    elif method == "joint":
        # Product-space distance combining signed and topological
        # d_diff. Two nodes are close iff they are close in both the
        # signed-Laplacian embedding (sign domains) and the unsigned-
        # Laplacian embedding (connectivity communities).
        #
        # d_joint(i,j)^2 = d_signed(i,j)^2 + d_topo(i,j)^2
        #
        # This is the Euclidean distance in the concatenated embedding
        # [V_signed * sqrt(w) | V_topo * sqrt(w)], which is the most
        # principled way to use both information sources without
        # hiding either.
        if G is None:
            raise ValueError("method='joint' requires the G parameter")
        # Signed channel
        D_signed = compute_signed_diffusion_distance(
            eigenvalues, eigenvectors_col, tau_star, use_abs=use_abs
        )
        # Topological channel: recompute eigendecomposition on |w|
        nodes = sorted(G.nodes())
        A_u = nx.adjacency_matrix(G, nodelist=nodes, weight="weight")
        A_u = np.abs(A_u.toarray().astype(float))
        D_u = np.diag(A_u.sum(axis=1))
        L_u = D_u - A_u
        eigv_u, eigV_u = np.linalg.eigh(L_u)
        D_topo = compute_signed_diffusion_distance(
            eigv_u, eigV_u, tau_star, use_abs=False
        )
        # Raw combination (no normalization). The relative magnitudes
        # of the two channels naturally reflect their information
        # content: on homogeneous lattices the topological channel is
        # small (nearly uniform), so the signed channel dominates —
        # correctly. On block-structured SBMs, both channels contribute.
        D = np.sqrt(D_signed ** 2 + D_topo ** 2)
        np.fill_diagonal(D, 0.0)
        D_cond = squareform(D)
        Z = linkage(D_cond, method=linkage_method)

        if n_clusters is not None:
            labels = fcluster(Z, t=n_clusters, criterion="maxclust")
        else:
            labels = _auto_k_silhouette(D, Z, N)

    elif method == "quantum_distance":
        D = compute_quantum_distance_matrix(eigenvectors_col, method="overlap")
        np.fill_diagonal(D, 0.0)
        D = np.maximum(D, D.T)
        mask = ~np.isfinite(D)
        if mask.any():
            D[mask] = np.nanmax(D[np.isfinite(D)]) * 2
        D_cond = squareform(D)
        Z = linkage(D_cond, method=linkage_method)

        if n_clusters is not None:
            labels = fcluster(Z, t=n_clusters, criterion="maxclust")
        else:
            labels = _auto_k_silhouette(D, Z, N)

    elif method == "eigenmode_sign":
        # Embed nodes using ALL eigenvectors, weighted by Boltzmann
        # factors at scale tau. Constant eigenvectors contribute zero
        # to pairwise distances automatically (v_i - v_j = 0), so no
        # mode skipping is needed — every mode carries geometry or
        # frustration information.
        #
        # On signed graphs, even the smallest-eigenvalue mode can be
        # highly informative (e.g. checkerboard on all-negative lattice).
        lam = np.abs(eigenvalues) if use_abs else eigenvalues
        boltz = np.exp(-tau_star * lam)
        w_modes = boltz / max(boltz.sum(), 1e-300)
        # Signed, weighted embedding over ALL modes
        V_emb = eigenvectors_col * np.sqrt(w_modes)[np.newaxis, :]
        from scipy.spatial.distance import pdist
        D_cond = pdist(V_emb, metric="euclidean")
        D = squareform(D_cond)
        Z = linkage(D_cond, method=linkage_method)

        if n_clusters is not None:
            labels = fcluster(Z, t=n_clusters, criterion="maxclust")
        else:
            labels = _auto_k_silhouette(D, Z, N)

    elif method == "ldos_entropy":
        S_local, _ = compute_ldos_entropy(
            eigenvalues, eigenvectors_col, tau_star, use_abs=use_abs
        )

        if n_clusters is None:
            from sklearn.cluster import KMeans
            from sklearn.metrics import silhouette_score

            best_k, best_score = 2, -1
            k_max = min(20, int(np.sqrt(N)))
            S_2d = S_local.reshape(-1, 1)
            for k in range(2, k_max + 1):
                km = KMeans(n_clusters=k, n_init=5, random_state=42)
                lab = km.fit_predict(S_2d)
                if len(np.unique(lab)) < 2:
                    continue
                score = silhouette_score(S_2d, lab)
                if score > best_score:
                    best_score = score
                    best_k = k
            n_clusters = best_k

        from sklearn.cluster import KMeans

        km = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)
        labels = km.fit_predict(S_local.reshape(-1, 1)) + 1
        D = None
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
        "linkage_matrix": Z,
        "distance_matrix": D,
    }


def _auto_k_silhouette(
    D: NDArray, Z: NDArray, N: int, k_max: int = 20
) -> NDArray:
    """Auto-detect cluster count via silhouette score on precomputed distance.

    Scans K = 2..k_max, picks the K with highest silhouette score,
    and returns the corresponding fcluster labels.
    """
    from sklearn.metrics import silhouette_score

    best_k, best_score = 2, -1.0
    k_upper = min(k_max, N // 2)
    for k in range(2, k_upper + 1):
        lab = fcluster(Z, t=k, criterion="maxclust")
        n_unique = len(np.unique(lab))
        if n_unique < 2 or n_unique >= N:
            continue
        score = silhouette_score(D, lab, metric="precomputed")
        if score > best_score:
            best_score = score
            best_k = k
    return fcluster(Z, t=best_k, criterion="maxclust")


def build_reduced_graph(
    G_original: nx.Graph,
    partition: list[set[int]],
    sign_method: str = "separate",
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
        How to assign reduced edge weights:

        - ``'separate'`` (default): Creates two edges between supernodes
          when both positive and negative cross-edges exist — one with
          the total positive weight, one with the total negative weight,
          both normalized by the number of pairs. Preserves frustration
          fraction faithfully.
        - ``'majority'``: Net signed weight / n_pairs. Simple but
          cancels mixed-sign edges, destroying frustration information.
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

    # Accumulate positive and negative cross-edges separately
    cross_pos = np.zeros((n_clusters, n_clusters))
    cross_neg = np.zeros((n_clusters, n_clusters))

    for u, v, d in G_original.edges(data=True):
        w = d.get("weight", 1.0)
        cu = node_to_cluster.get(u)
        cv = node_to_cluster.get(v)
        if cu is not None and cv is not None:
            if w >= 0:
                cross_pos[cu, cv] += w
                cross_pos[cv, cu] += w
            else:
                cross_neg[cu, cv] += w
                cross_neg[cv, cu] += w

    cross_weights = cross_pos + cross_neg

    # Build reduced graph
    G_reduced = nx.Graph()
    for c_idx in range(n_clusters):
        G_reduced.add_node(
            c_idx,
            size=len(partition[c_idx]),
            internal_weight=cross_weights[c_idx, c_idx] / 2,
            original_nodes=partition[c_idx],
        )

    for i in range(n_clusters):
        for j in range(i + 1, n_clusters):
            n_pairs = len(partition[i]) * len(partition[j])
            if n_pairs == 0:
                continue

            w_pos = cross_pos[i, j]
            w_neg = cross_neg[i, j]  # negative number

            if sign_method == "separate":
                # Dominant sign determines reduced edge weight;
                # magnitude proportional to edge density.
                # Frustration info stored as node/edge attributes.
                w_total = w_pos + w_neg  # net signed
                n_edges = (w_pos + abs(w_neg))  # total |weight|
                if n_edges < 1e-10:
                    continue
                frac_neg = abs(w_neg) / n_edges
                w_norm = n_edges / n_pairs  # density of connection
                if w_norm > min_edge_weight:
                    sign = 1.0 if w_pos >= abs(w_neg) else -1.0
                    G_reduced.add_edge(
                        i, j,
                        weight=sign * w_norm,
                        frustration=frac_neg,
                        n_pos=int(w_pos),
                        n_neg=int(abs(w_neg)),
                    )
            elif sign_method == "majority":
                w_raw = w_pos + w_neg
                w_norm = abs(w_raw) / n_pairs
                if w_norm > min_edge_weight and abs(w_raw) > 1e-10:
                    G_reduced.add_edge(i, j, weight=w_raw / n_pairs)
            else:
                raise ValueError(f"Unknown sign_method: {sign_method}")

    # Metadata
    metadata = {
        "node_map": {i: partition[i] for i in range(n_clusters)},
        "cluster_sizes": [len(p) for p in partition],
        "internal_weights": [cross_weights[i, i] / 2 for i in range(n_clusters)],
        "cross_weight_matrix": cross_weights,
        "cross_pos": cross_pos,
        "cross_neg": cross_neg,
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
    partition_method: str = "joint",
    sign_method: str = "separate",
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
        Partition method. Default ``'eigenmode_sign'`` (clusters by
        low-mode sign pattern — groups nodes into frustration-consistent
        domains, keeping frustration at cluster boundaries as an RG
        flow variable). Use ``'topological'`` for block-structure
        preserving coarse-graining.
    sign_method : str
        Edge construction: ``'separate'`` (default, preserves per-edge
        frustration info) or ``'majority'`` (net signed weight).
    n_clusters : int, optional
        Force cluster count. If None, auto-detect via silhouette scan.
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
        G=G,
    )

    # Build reduced graph
    G_reduced, metadata = build_reduced_graph(
        G, part_info["partition"], sign_method=sign_method
    )

    # Compute intra-cluster frustration
    intra_frust = []
    for cluster in part_info["partition"]:
        sub = G.subgraph(list(cluster))
        ws = [d.get("weight", 1) for _, _, d in sub.edges(data=True)]
        if ws:
            intra_frust.append(sum(1 for w in ws if w < 0) / len(ws))
        else:
            intra_frust.append(0.0)
    metadata["intra_cluster_frustration"] = np.array(intra_frust)
    metadata["mean_intra_frustration"] = float(np.mean(intra_frust))

    # Compute boundary frustration (per-edge frustration of reduced edges)
    edge_frusts = [
        d.get("frustration", 1.0 if d.get("weight", 1) < 0 else 0.0)
        for _, _, d in G_reduced.edges(data=True)
    ]
    if edge_frusts:
        metadata["mean_boundary_frustration"] = float(np.mean(edge_frusts))
    else:
        metadata["mean_boundary_frustration"] = 0.0

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
    partition_method: str = "joint",
    sign_method: str = "separate",
    use_abs: bool = True,
    verbose: bool = False,
) -> list[dict[str, Any]]:
    """Iterate spectral RG steps to produce a coarse-graining flow.

    Frustration is tracked as an RG flow variable at each step:
    it concentrates at cluster boundaries under coarse-graining
    (the eigenmode partition creates internally consistent domains).
    Below p_c frustration flows toward zero; above p_c it persists
    and may increase — reflecting the irreducible disorder.

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
        Default ``'eigenmode_sign'`` (frustration-preserving).
        ``'topological'`` for block-structure coarse-graining.
    sign_method : str
        Edge construction: ``'separate'`` (default) or ``'majority'``.
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
            sign_method=sign_method,
            use_abs=use_abs,
        )

        flow.append(step)

        if verbose:
            meta = step["metadata"]
            ff = meta["frustration_fraction"]
            bf = meta.get("mean_boundary_frustration", 0)
            inf = meta.get("mean_intra_frustration", 0)
            print(
                f"  Step {step_count}: {step['N_original']} → "
                f"{step['N_reduced']} nodes | τ*={step['tau_star']:.3f} "
                f"| f_edge={ff:.3f} f_boundary={bf:.3f} f_intra={inf:.3f}"
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
        Arrays indexed by RG step:

        - ``n_nodes``: node count at each scale (length = n_steps + 1)
        - ``n_edges``: edge count of reduced graph
        - ``frustration_fraction``: fraction of negative edges in reduced graph
        - ``boundary_frustration``: mean per-edge frustration at cluster
          boundaries — the RG flow variable for frustration
        - ``intra_frustration``: mean intra-cluster frustration (should
          stay low if eigenmode partition is working)
        - ``spectral_gap``: smallest nonzero |eigenvalue|
        - ``tau_stars``: RG scale at each step
        - ``entropy_at_tau_star``: global entropy at the chosen scale
    """
    n_nodes = []
    n_edges = []
    frustration = []
    boundary_frustration = []
    intra_frustration = []
    spectral_gap = []
    tau_stars = []
    entropy_at_tau_star = []

    for step in flow:
        G_r = step["G_reduced"]
        eigv_r = step["eigenvalues_reduced"]
        meta = step["metadata"]

        n_nodes.append(step["N_original"])
        n_edges.append(G_r.number_of_edges())
        frustration.append(meta["frustration_fraction"])
        boundary_frustration.append(
            meta.get("mean_boundary_frustration", 0.0)
        )
        intra_frustration.append(
            meta.get("mean_intra_frustration", 0.0)
        )

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
        "boundary_frustration": np.array(boundary_frustration),
        "intra_frustration": np.array(intra_frustration),
        "spectral_gap": np.array(spectral_gap),
        "tau_stars": np.array(tau_stars),
        "entropy_at_tau_star": np.array(entropy_at_tau_star),
    }
