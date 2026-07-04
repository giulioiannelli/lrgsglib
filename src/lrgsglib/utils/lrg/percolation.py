"""Percolation analysis for signed planar graphs.

Two complementary methods for extracting percolation thresholds
(surface-code error thresholds) from signed graphs:

**Eigenmode sign-cluster method** (arXiv:2504.00144):
Binarize the lowest Laplacian eigenvector by sign, extract connected
components, measure giant-cluster fraction P_inf vs disorder p.

**Direct defect identification** (Dennis et al. 2002):
Enumerate frustrated plaquettes (Z-syndromes) and frustrated vertices
(X-syndromes), build syndrome graphs, measure their percolation.
Geometry-agnostic: works on any planar graph without constructing
the dual.

In both cases, the susceptibility peak ``chi = sqrt(N) * std(P_inf)``
locates the threshold p_c.
"""

from __future__ import annotations

from typing import Any, Callable, Optional, TYPE_CHECKING

import networkx as nx
import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from ...graphs.nx.SignedGraphNX.SignedGraphNX import SignedGraphNX

__all__ = [
    # --- eigenmode method ---
    "percolation_sweep",
    "find_threshold",
    "dual_percolation_sweep",
    "compute_pareto_point",
    # --- generalized dual (any graph) ---
    "build_cycle_dual",
    "build_signed_cycle_dual",
    "dual_percolation_sweep_general",
    "compute_pareto_point_general",
    # --- direct defect method ---
    "identify_frustrated_faces",
    "identify_frustrated_vertices",
    "giant_component_fraction",
    "build_face_syndrome",
    "build_vertex_syndrome",
    "defect_percolation_sweep",
    "compute_pareto_point_direct",
]


def percolation_sweep(
    graph_factory: Callable[..., Any],
    graph_params: dict[str, Any],
    p_range: NDArray[np.floating],
    n_realizations: int = 50,
    which: int = 0,
    seed_base: int = 0,
    backend: str = "scipy",
    verbose: bool = False,
) -> dict[str, Any]:
    """Sweep disorder strength, measuring giant-cluster fraction.

    For each *p* in ``p_range``, creates ``n_realizations`` independent
    disorder realizations, computes the lowest eigenvector, binarizes it
    by sign, and measures the giant-cluster fraction P_inf.

    Parameters
    ----------
    graph_factory : callable
        Class or function returning a ``SignedGraphNX``.  Must accept
        ``pflip`` and ``seed`` keyword arguments (e.g. ``Lattice2DNX``).
    graph_params : dict
        Fixed construction parameters (e.g. ``{'side1': 32, 'geo': 'sqr',
        'pbc': False}``).
    p_range : ndarray
        Array of ``pflip`` values to sweep.
    n_realizations : int
        Disorder realizations per p value.
    which : int
        Eigenvector index (0 = ground state).
    seed_base : int
        Base seed for reproducibility.
    backend : str
        Eigendecomposition backend (``'scipy'`` recommended for sparse).
    verbose : bool
        Print progress.

    Returns
    -------
    dict
        ``'p'``: p_range,
        ``'Pinf_mean'``: mean P_inf(p),
        ``'Pinf_std'``: std P_inf(p),
        ``'chi'``: susceptibility sqrt(N) * std(P_inf),
        ``'N'``: system size,
        ``'n_realizations'``: realizations count,
        ``'Pinf_all'``: full (n_p, n_real) array.
    """
    n_p = len(p_range)
    Pinf_all = np.zeros((n_p, n_realizations))
    N = None

    for ip, p in enumerate(p_range):
        if verbose:
            print(f"  p={p:.4f} ({ip+1}/{n_p})", end="\r", flush=True)
        for ir in range(n_realizations):
            seed = seed_base + ip * n_realizations + ir
            sg = graph_factory(**graph_params, pflip=p, seed=seed)
            sg.flip_random_fract_edges()
            sg.compute_k_eigvV(k=which + 1, backend=backend)
            sg.compute_pinf(which=which)
            Pinf_all[ip, ir] = sg.Pinf
            if N is None:
                N = sg.N

    Pinf_mean = Pinf_all.mean(axis=1)
    Pinf_std = Pinf_all.std(axis=1)
    chi = np.sqrt(N) * Pinf_std

    if verbose:
        print()

    return {
        "p": np.asarray(p_range),
        "Pinf_mean": Pinf_mean,
        "Pinf_std": Pinf_std,
        "chi": chi,
        "N": N,
        "n_realizations": n_realizations,
        "Pinf_all": Pinf_all,
    }


def find_threshold(
    p_range: NDArray[np.floating],
    chi_values: NDArray[np.floating],
    method: str = "peak",
) -> float:
    """Locate the percolation threshold from susceptibility data.

    Parameters
    ----------
    p_range : ndarray
        Disorder-strength values.
    chi_values : ndarray
        Susceptibility values (same length as *p_range*).
    method : str
        ``'peak'``: raw argmax.
        ``'fit'``: quadratic fit near the peak for sub-grid resolution.

    Returns
    -------
    float
        Estimated p_c.
    """
    if method == "peak":
        return float(p_range[np.argmax(chi_values)])

    if method == "fit":
        idx_peak = int(np.argmax(chi_values))
        lo = max(0, idx_peak - 2)
        hi = min(len(p_range), idx_peak + 3)
        p_local = p_range[lo:hi]
        chi_local = chi_values[lo:hi]
        if len(p_local) < 3:
            return float(p_range[idx_peak])
        coeffs = np.polyfit(p_local, chi_local, 2)
        # Vertex of parabola: -b / (2a)
        if coeffs[0] != 0:
            p_c = -coeffs[1] / (2 * coeffs[0])
            if p_local[0] <= p_c <= p_local[-1]:
                return float(p_c)
        return float(p_range[idx_peak])

    raise ValueError(f"Unknown method: {method!r}")


def dual_percolation_sweep(
    graph_factory: Callable[..., Any],
    graph_params: dict[str, Any],
    p_range: NDArray[np.floating],
    n_realizations: int = 50,
    which: int = 0,
    seed_base: int = 0,
    backend: str = "scipy",
    include_outer_face: bool = False,
    verbose: bool = False,
) -> dict[str, dict[str, Any]]:
    """Percolation sweep on both a graph and its planar dual.

    For each realization builds the primal graph with disorder, then
    constructs the planar dual with inherited signs, and measures
    P_inf on both.

    Parameters
    ----------
    graph_factory : callable
        Class returning a ``SignedGraphNX`` (must be planar when
        ``pbc=False``).
    graph_params : dict
        Must include ``pbc=False`` for planarity.
    p_range : ndarray
        Disorder-strength values.
    n_realizations : int
        Realizations per p value.
    which : int
        Eigenvector index.
    seed_base : int
        Base seed.
    backend : str
        Eigen-backend.
    include_outer_face : bool
        Whether the dual includes the infinite face.
    verbose : bool
        Print progress.

    Returns
    -------
    dict
        ``'primal'``: percolation_sweep result dict,
        ``'dual'``: matching result dict for the planar dual.
    """
    from ...graphs.nx.funcs.dual import build_signed_dual

    n_p = len(p_range)
    N_primal = None
    N_dual = None
    Pinf_primal = np.zeros((n_p, n_realizations))
    Pinf_dual = np.zeros((n_p, n_realizations))

    for ip, p in enumerate(p_range):
        if verbose:
            print(f"  p={p:.4f} ({ip+1}/{n_p})", end="\r", flush=True)
        for ir in range(n_realizations):
            seed = seed_base + ip * n_realizations + ir
            sg = graph_factory(**graph_params, pflip=p, seed=seed)
            sg.flip_random_fract_edges()

            # Primal
            sg.compute_k_eigvV(k=which + 1, backend=backend)
            sg.compute_pinf(which=which)
            Pinf_primal[ip, ir] = sg.Pinf

            # Dual
            dual_sg = build_signed_dual(
                sg, include_outer_face=include_outer_face
            )
            dual_sg.compute_k_eigvV(k=which + 1, backend=backend)
            dual_sg.compute_pinf(which=which)
            Pinf_dual[ip, ir] = dual_sg.Pinf

            if N_primal is None:
                N_primal = sg.N
                N_dual = dual_sg.N

    if verbose:
        print()

    def _pack(Pinf_arr: NDArray, N: int) -> dict[str, Any]:
        mean = Pinf_arr.mean(axis=1)
        std = Pinf_arr.std(axis=1)
        return {
            "p": np.asarray(p_range),
            "Pinf_mean": mean,
            "Pinf_std": std,
            "chi": np.sqrt(N) * std,
            "N": N,
            "n_realizations": n_realizations,
            "Pinf_all": Pinf_arr,
        }

    return {
        "primal": _pack(Pinf_primal, N_primal),
        "dual": _pack(Pinf_dual, N_dual),
    }


def compute_pareto_point(
    graph_factory: Callable[..., Any],
    graph_params: dict[str, Any],
    p_range: NDArray[np.floating],
    n_realizations: int = 50,
    which: int = 0,
    seed_base: int = 0,
    backend: str = "scipy",
    include_outer_face: bool = False,
    threshold_method: str = "fit",
    verbose: bool = False,
) -> tuple[float, float, dict[str, dict[str, Any]]]:
    """Compute the (p_X, p_Z) Pareto point for a lattice geometry.

    Runs :func:`dual_percolation_sweep` and extracts thresholds from
    both primal and dual susceptibility peaks.

    Parameters
    ----------
    graph_factory, graph_params, p_range, n_realizations, which, seed_base, backend, include_outer_face, verbose :
        Forwarded to :func:`dual_percolation_sweep`.
    threshold_method : str
        Forwarded to :func:`find_threshold`.

    Returns
    -------
    p_X : float
        X-error threshold (from dual susceptibility peak).
    p_Z : float
        Z-error threshold (from primal susceptibility peak).
    sweep_data : dict
        Full sweep data (``'primal'`` and ``'dual'`` dicts).
    """
    sweep = dual_percolation_sweep(
        graph_factory=graph_factory,
        graph_params=graph_params,
        p_range=p_range,
        n_realizations=n_realizations,
        which=which,
        seed_base=seed_base,
        backend=backend,
        include_outer_face=include_outer_face,
        verbose=verbose,
    )

    p_Z = find_threshold(
        sweep["primal"]["p"], sweep["primal"]["chi"], method=threshold_method
    )
    p_X = find_threshold(
        sweep["dual"]["p"], sweep["dual"]["chi"], method=threshold_method
    )

    return p_X, p_Z, sweep


# ---------------------------------------------------------------
# Direct defect identification (geometry-agnostic)
# ---------------------------------------------------------------


def identify_frustrated_faces(
    G: "nx.Graph",
    face_info: dict[str, Any],
    weight: str = "weight",
) -> list[int]:
    """Identify frustrated plaquettes (Z-defect syndromes).

    A face is *frustrated* when it has an odd number of negative
    boundary edges.  This is the Z-syndrome of the surface code:
    a frustrated plaquette flags a Z-error chain crossing it.

    Parameters
    ----------
    G : nx.Graph
        Primal graph with signed edge weights.
    face_info : dict
        Output of :func:`~lrgsglib.graphs.nx.funcs.dual.build_planar_dual`,
        containing ``'faces'`` (list of node-tuples).
    weight : str
        Edge attribute storing the sign.

    Returns
    -------
    list[int]
        Indices of frustrated faces.
    """
    frustrated: list[int] = []
    for fi, face in enumerate(face_info["faces"]):
        n = len(face)
        n_neg = 0
        for k in range(n):
            u, v = face[k], face[(k + 1) % n]
            if G.has_edge(u, v):
                w = G[u][v].get(weight, 1)
                if w < 0:
                    n_neg += 1
        if n_neg % 2 == 1:
            frustrated.append(fi)
    return frustrated


def identify_frustrated_vertices(
    G: "nx.Graph",
    weight: str = "weight",
) -> list:
    """Identify frustrated vertices (X-defect syndromes).

    A vertex is *frustrated* when it has an odd number of negative
    incident edges.  This is the X-syndrome of the surface code:
    a frustrated vertex flags an X-error chain crossing it.

    Parameters
    ----------
    G : nx.Graph
        Graph with signed edge weights.
    weight : str
        Edge attribute storing the sign.

    Returns
    -------
    list
        Node labels of frustrated vertices.
    """
    frustrated = []
    for v in G.nodes():
        n_neg = sum(
            1 for u in G.neighbors(v) if G[v][u].get(weight, 1) < 0
        )
        if n_neg % 2 == 1:
            frustrated.append(v)
    return frustrated


def giant_component_fraction(
    syndrome: "nx.Graph",
) -> float:
    """Giant-component fraction of a syndrome graph.

    Parameters
    ----------
    syndrome : nx.Graph
        Syndrome graph (nodes = defects, edges = adjacency).

    Returns
    -------
    float
        Fraction of nodes in the largest connected component,
        or 0.0 if the graph is empty.
    """
    n = syndrome.number_of_nodes()
    if n == 0:
        return 0.0
    largest = max(len(c) for c in nx.connected_components(syndrome))
    return largest / n


def build_face_syndrome(
    face_info: dict[str, Any],
    frustrated_faces: list[int],
) -> nx.Graph:
    """Build syndrome graph connecting adjacent frustrated faces.

    Two frustrated faces are *adjacent* if they share a primal edge.
    The resulting graph has one node per frustrated face and edges
    between face-pairs that share an edge.

    Parameters
    ----------
    face_info : dict
        Output of :func:`~lrgsglib.graphs.nx.funcs.dual.build_planar_dual`.
    frustrated_faces : list[int]
        Face indices identified by :func:`identify_frustrated_faces`.

    Returns
    -------
    nx.Graph
        Syndrome graph (nodes = frustrated face indices).
    """
    fset = set(frustrated_faces)
    syn = nx.Graph()
    syn.add_nodes_from(frustrated_faces)

    for _e, fp in face_info["edge_to_faces"].items():
        if len(fp) == 2:
            f1, f2 = fp
            if f1 in fset and f2 in fset:
                syn.add_edge(f1, f2)
    return syn


def build_vertex_syndrome(
    G: nx.Graph,
    frustrated_vertices: list,
) -> nx.Graph:
    """Build syndrome graph: subgraph of *G* induced by frustrated vertices.

    Two frustrated vertices are connected if they share an edge in *G*.

    Parameters
    ----------
    G : nx.Graph
        The primal graph.
    frustrated_vertices : list
        Vertex labels from :func:`identify_frustrated_vertices`.

    Returns
    -------
    nx.Graph
        Syndrome graph (nodes = frustrated vertex labels).
    """
    return G.subgraph(frustrated_vertices).copy()


def defect_percolation_sweep(
    graph_factory: Callable[..., Any],
    graph_params: dict[str, Any],
    p_range: NDArray[np.floating],
    n_realizations: int = 50,
    seed_base: int = 0,
    include_outer_face: bool = False,
    verbose: bool = False,
) -> dict[str, dict[str, Any]]:
    """Sweep disorder, measuring percolation of Z- and X-defect syndromes.

    For each *p*, creates ``n_realizations`` disorder instances,
    enumerates frustrated plaquettes (Z-syndromes) and frustrated
    vertices (X-syndromes), builds syndrome graphs, and measures their
    giant-component fractions.

    Works on **any** planar graph (regular lattices, fractals, random
    planar graphs) without constructing the dual.

    .. note::

       Syndrome percolation measures a different observable than the
       RBIM phase transition.  For actual error thresholds (p_X, p_Z),
       use :func:`compute_pareto_point` (eigenmode method).  This
       function provides **diagnostic analysis** of the defect
       landscape: defect densities, syndrome clustering, and
       face-size effects.

    Parameters
    ----------
    graph_factory : callable
        Class returning a ``SignedGraphNX`` (must be planar).
    graph_params : dict
        Fixed construction parameters.
    p_range : ndarray
        Disorder strengths to sweep.
    n_realizations : int
        Independent disorder samples per p.
    seed_base : int
        Base seed.
    include_outer_face : bool
        Whether to include the outer face in Z-syndrome identification.
    verbose : bool
        Print progress.

    Returns
    -------
    dict
        ``'Z'``: dict with keys ``'p'``, ``'Pinf_mean'``, ``'Pinf_std'``,
        ``'chi'``, ``'N'``, ``'n_defects_mean'``, ``'Pinf_all'``.
        ``'X'``: dict with same keys.
        ``'face_info_sample'``: face_info from first realization (for
        reference / face-size analysis).
    """
    from ...graphs.nx.funcs.dual import build_planar_dual

    n_p = len(p_range)
    Pinf_Z = np.zeros((n_p, n_realizations))
    Pinf_X = np.zeros((n_p, n_realizations))
    n_defects_Z = np.zeros((n_p, n_realizations))
    n_defects_X = np.zeros((n_p, n_realizations))
    N_vertices: Optional[int] = None
    N_faces: Optional[int] = None
    face_info_sample: Optional[dict] = None

    for ip, p in enumerate(p_range):
        if verbose:
            print(f"  p={p:.4f} ({ip+1}/{n_p})", end="\r", flush=True)
        for ir in range(n_realizations):
            seed = seed_base + ip * n_realizations + ir
            sg = graph_factory(**graph_params, pflip=p, seed=seed)
            sg.flip_random_fract_edges()
            G = sg.gr[sg.on_g]

            # Face enumeration via planar embedding
            _, fi = build_planar_dual(G, include_outer_face=True)

            # Z-syndromes: frustrated plaquettes
            frust_faces = identify_frustrated_faces(G, fi)
            if not include_outer_face:
                frust_faces = [
                    f for f in frust_faces
                    if f != fi["outer_face_idx"]
                ]

            syn_Z = build_face_syndrome(fi, frust_faces)
            Pinf_Z[ip, ir] = giant_component_fraction(syn_Z)
            n_defects_Z[ip, ir] = len(frust_faces)

            # X-syndromes: frustrated vertices
            frust_verts = identify_frustrated_vertices(G)
            syn_X = build_vertex_syndrome(G, frust_verts)
            Pinf_X[ip, ir] = giant_component_fraction(syn_X)
            n_defects_X[ip, ir] = len(frust_verts)

            if N_vertices is None:
                N_vertices = G.number_of_nodes()
                n_total_faces = len(fi["faces"])
                if not include_outer_face:
                    n_total_faces -= 1
                N_faces = n_total_faces
                face_info_sample = fi

    if verbose:
        print()

    def _pack(
        Pinf_arr: NDArray, N: int, n_def: NDArray
    ) -> dict[str, Any]:
        mean = Pinf_arr.mean(axis=1)
        std = Pinf_arr.std(axis=1)
        return {
            "p": np.asarray(p_range),
            "Pinf_mean": mean,
            "Pinf_std": std,
            "chi": np.sqrt(N) * std,
            "N": N,
            "n_realizations": n_realizations,
            "n_defects_mean": n_def.mean(axis=1),
            "Pinf_all": Pinf_arr,
        }

    result: dict[str, Any] = {
        "Z": _pack(Pinf_Z, N_faces, n_defects_Z),
        "X": _pack(Pinf_X, N_vertices, n_defects_X),
        "face_info_sample": face_info_sample,
    }
    return result


def compute_pareto_point_direct(
    graph_factory: Callable[..., Any],
    graph_params: dict[str, Any],
    p_range: NDArray[np.floating],
    n_realizations: int = 50,
    seed_base: int = 0,
    include_outer_face: bool = False,
    threshold_method: str = "fit",
    verbose: bool = False,
) -> tuple[float, float, dict[str, Any]]:
    """Compute (p_X, p_Z) via direct defect percolation.

    Runs :func:`defect_percolation_sweep` and extracts thresholds
    from both Z-syndrome and X-syndrome susceptibility peaks.

    Returns
    -------
    p_X : float
        X-error threshold (frustrated-vertex syndrome percolation).
    p_Z : float
        Z-error threshold (frustrated-plaquette syndrome percolation).
    sweep_data : dict
        Full sweep data with ``'Z'`` and ``'X'`` dicts.
    """
    sweep = defect_percolation_sweep(
        graph_factory=graph_factory,
        graph_params=graph_params,
        p_range=p_range,
        n_realizations=n_realizations,
        seed_base=seed_base,
        include_outer_face=include_outer_face,
        verbose=verbose,
    )

    p_Z = find_threshold(
        sweep["Z"]["p"], sweep["Z"]["chi"], method=threshold_method
    )
    p_X = find_threshold(
        sweep["X"]["p"], sweep["X"]["chi"], method=threshold_method
    )

    return p_X, p_Z, sweep


def default_p_range(
    p_min: float = 0.005,
    p_max: float = 0.5,
    n_fine: int = 25,
    n_coarse: int = 8,
    focus: Optional[tuple[float, float]] = None,
) -> NDArray[np.floating]:
    """Build a p-range with higher density near the expected transition.

    Parameters
    ----------
    p_min, p_max : float
        Sweep bounds.
    n_fine : int
        Points in the fine-grained region.
    n_coarse : int
        Points in each coarse tail.
    focus : tuple[float, float], optional
        ``(lo, hi)`` bounds for the fine-grained region.
        Defaults to ``(0.06, 0.25)``.

    Returns
    -------
    ndarray
        Sorted unique p values.
    """
    if focus is None:
        focus = (0.06, 0.25)
    p1 = np.linspace(p_min, focus[0], n_coarse, endpoint=False)
    p2 = np.linspace(focus[0], focus[1], n_fine)
    p3 = np.linspace(focus[1], p_max, n_coarse + 1)[1:]
    return np.unique(np.concatenate([p1, p2, p3]))


# ---------------------------------------------------------------
# Generalized dual via minimum cycle basis
# ---------------------------------------------------------------


def build_cycle_dual(
    G: nx.Graph,
    weight: str = "weight",
) -> tuple[nx.Graph, dict[str, Any]]:
    """Build the generalized dual of any connected graph.

    Uses the **minimum cycle basis** (MCB) to identify independent
    cycles (= generalized plaquettes).  Two cycles are adjacent in the
    dual if they share at least one primal edge.

    For **planar** graphs this is equivalent to
    :func:`~lrgsglib.graphs.nx.funcs.dual.build_planar_dual`
    (verified: 4×4 grid gives 9 nodes / 12 edges in both).

    Parameters
    ----------
    G : nx.Graph
        A connected graph with signed edge weights.
    weight : str
        Edge attribute storing the sign.

    Returns
    -------
    dual : nx.Graph
        Nodes = cycle indices, edges between cycles sharing ≥1 primal
        edge.  Edge weights = product of shared primal-edge signs.
    cycle_info : dict
        ``'cycles'``: list of node-tuples (one per independent cycle),
        ``'edge_to_cycles'``: ``{(u,v): [cycle indices]}``,
        ``'cycle_sizes'``: list of cycle perimeters,
        ``'n_cycles'``: number of independent cycles (= b₁).
    """
    mcb = nx.minimum_cycle_basis(G)
    n_cycles = len(mcb)

    # Map edges → list of containing cycle indices
    edge_to_cycles: dict[tuple, list[int]] = {}
    for ci, cycle in enumerate(mcb):
        n = len(cycle)
        for k in range(n):
            u, v = cycle[k], cycle[(k + 1) % n]
            e = (min(u, v), max(u, v))
            edge_to_cycles.setdefault(e, []).append(ci)

    # Build dual graph
    dual = nx.Graph()
    dual.add_nodes_from(range(n_cycles))
    for ci in range(n_cycles):
        dual.nodes[ci]["cycle_size"] = len(mcb[ci])

    # Connect cycles that share edges; inherit sign
    adjacency_signs: dict[tuple[int, int], list[int]] = {}
    for e, clist in edge_to_cycles.items():
        if len(clist) < 2:
            continue
        u, v = e
        w = G[u][v].get(weight, 1)
        for i in range(len(clist)):
            for j in range(i + 1, len(clist)):
                c1, c2 = min(clist[i], clist[j]), max(clist[i], clist[j])
                adjacency_signs.setdefault((c1, c2), []).append(w)

    for (c1, c2), signs in adjacency_signs.items():
        # Product of shared edge signs
        prod = 1
        for s in signs:
            prod *= (1 if s >= 0 else -1)
        dual.add_edge(c1, c2, **{weight: prod})

    cycle_info = {
        "cycles": [tuple(c) for c in mcb],
        "edge_to_cycles": edge_to_cycles,
        "cycle_sizes": [len(c) for c in mcb],
        "n_cycles": n_cycles,
    }
    return dual, cycle_info


def build_signed_cycle_dual(
    sg: "SignedGraphNX",
    on_g: Optional[str] = None,
) -> "SignedGraphNX":
    """Build a :class:`SignedGraphNX` wrapping the cycle dual.

    Generalizes :func:`~lrgsglib.graphs.nx.funcs.dual.build_signed_dual`
    to arbitrary (non-planar) graphs.

    Parameters
    ----------
    sg : SignedGraphNX
        Source signed graph.
    on_g : str, optional
        Graph representation key.

    Returns
    -------
    SignedGraphNX
        Dual signed graph ready for ``compute_k_eigvV`` /
        ``compute_pinf``.
    """
    from ...graphs.nx.SignedGraphNX.SignedGraphNX import SignedGraphNX
    from ...config.const import SG_REPR

    if on_g is None:
        on_g = getattr(sg, "on_g", SG_REPR)

    G = sg.gr[on_g]
    dual_nx, cycle_info = build_cycle_dual(G)

    if dual_nx.number_of_nodes() == 0:
        # Degenerate: tree graph has no cycles
        dual_sg = SignedGraphNX(nx.Graph(), pflip=0.0, make_dir_tree=False)
        dual_sg._cycle_dual_info = cycle_info
        return dual_sg

    dual_nx = nx.convert_node_labels_to_integers(dual_nx)

    dual_sg = SignedGraphNX(
        dual_nx,
        pflip=0.0,
        make_dir_tree=False,
    )
    dual_sg.upd_edge_sets(on_g=dual_sg.on_g)
    dual_sg.pflip = (
        len(dual_sg.fleset[dual_sg.on_g]) / dual_sg.Ne
        if dual_sg.Ne > 0
        else 0.0
    )
    dual_sg.Ne_flips = len(dual_sg.fleset[dual_sg.on_g])
    dual_sg._cycle_dual_info = cycle_info
    dual_sg._cycle_dual_of = sg

    return dual_sg


def dual_percolation_sweep_general(
    graph_factory: Callable[..., Any],
    graph_params: dict[str, Any],
    p_range: NDArray[np.floating],
    n_realizations: int = 50,
    which: int = 0,
    seed_base: int = 0,
    backend: str = "scipy",
    verbose: bool = False,
) -> dict[str, dict[str, Any]]:
    """Percolation sweep on both a graph and its cycle dual.

    Works on **any** connected graph (planar or non-planar).
    Uses the minimum cycle basis to construct the generalized dual.

    Parameters
    ----------
    graph_factory : callable
        Class returning a ``SignedGraphNX``.
    graph_params : dict
        Construction parameters.
    p_range : ndarray
        Disorder-strength values.
    n_realizations : int
        Realizations per p value.
    which : int
        Eigenvector index.
    seed_base : int
        Base seed.
    backend : str
        Eigen-backend.
    verbose : bool
        Print progress.

    Returns
    -------
    dict
        ``'primal'``: sweep result dict (P_inf on G → p_Z),
        ``'dual'``: sweep result dict (P_inf on cycle dual → p_X).
    """
    n_p = len(p_range)
    N_primal: Optional[int] = None
    N_dual: Optional[int] = None
    Pinf_primal = np.zeros((n_p, n_realizations))
    Pinf_dual = np.zeros((n_p, n_realizations))

    for ip, p in enumerate(p_range):
        if verbose:
            print(f"  p={p:.4f} ({ip+1}/{n_p})", end="\r", flush=True)
        for ir in range(n_realizations):
            seed = seed_base + ip * n_realizations + ir
            sg = graph_factory(**graph_params, pflip=p, seed=seed)
            sg.flip_random_fract_edges()

            # Primal
            sg.compute_k_eigvV(k=which + 1, backend=backend)
            sg.compute_pinf(which=which)
            Pinf_primal[ip, ir] = sg.Pinf

            # Cycle dual
            dual_sg = build_signed_cycle_dual(sg)
            if dual_sg.N > 1:
                dual_sg.compute_k_eigvV(k=which + 1, backend=backend)
                dual_sg.compute_pinf(which=which)
                Pinf_dual[ip, ir] = dual_sg.Pinf
            else:
                Pinf_dual[ip, ir] = 0.0

            if N_primal is None:
                N_primal = sg.N
                N_dual = dual_sg.N

    if verbose:
        print()

    def _pack(Pinf_arr: NDArray, N: int) -> dict[str, Any]:
        mean = Pinf_arr.mean(axis=1)
        std = Pinf_arr.std(axis=1)
        return {
            "p": np.asarray(p_range),
            "Pinf_mean": mean,
            "Pinf_std": std,
            "chi": np.sqrt(max(N, 1)) * std,
            "N": N,
            "n_realizations": n_realizations,
            "Pinf_all": Pinf_arr,
        }

    return {
        "primal": _pack(Pinf_primal, N_primal),
        "dual": _pack(Pinf_dual, N_dual),
    }


def compute_pareto_point_general(
    graph_factory: Callable[..., Any],
    graph_params: dict[str, Any],
    p_range: NDArray[np.floating],
    n_realizations: int = 50,
    which: int = 0,
    seed_base: int = 0,
    backend: str = "scipy",
    threshold_method: str = "fit",
    verbose: bool = False,
) -> tuple[float, float, dict[str, dict[str, Any]]]:
    """Compute (p_X, p_Z) for any graph via cycle dual.

    Generalizes :func:`compute_pareto_point` to non-planar graphs.

    Returns
    -------
    p_X : float
        X-error threshold (cycle dual susceptibility peak).
    p_Z : float
        Z-error threshold (primal susceptibility peak).
    sweep_data : dict
        Full sweep data.
    """
    sweep = dual_percolation_sweep_general(
        graph_factory=graph_factory,
        graph_params=graph_params,
        p_range=p_range,
        n_realizations=n_realizations,
        which=which,
        seed_base=seed_base,
        backend=backend,
        verbose=verbose,
    )

    p_Z = find_threshold(
        sweep["primal"]["p"], sweep["primal"]["chi"],
        method=threshold_method,
    )
    p_X = find_threshold(
        sweep["dual"]["p"], sweep["dual"]["chi"],
        method=threshold_method,
    )

    return p_X, p_Z, sweep
