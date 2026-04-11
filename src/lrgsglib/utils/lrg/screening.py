"""Topology screening for quantum error correction.

Systematic infrastructure for evaluating planar graph topologies
as surface-code substrates.  Each topology is scored by:

1. **Face uniformity** — coefficient of variation of face sizes.
   Uniform faces (CV ≈ 0) are required for viable Z-stabilizers.
2. **Error thresholds** — (p_X, p_Z) from the eigenmode method
   (primal + dual susceptibility peaks).
3. **Distance above the Pareto front** — signed distance from the
   Fujii & Tokunaga (2012) regular-lattice tradeoff curve.

The screening pipeline uses a two-stage filter: cheap planarity +
face-uniformity check first, expensive eigenmode sweep only on
surviving candidates.
"""

from __future__ import annotations

import time
from collections import Counter
from typing import Any, Callable, Optional

import networkx as nx
import numpy as np
import pandas as pd
from numpy.typing import NDArray

from .percolation import (
    compute_pareto_point,
    compute_pareto_point_general,
    default_p_range,
    find_threshold,
    percolation_sweep,
)

__all__ = [
    "compute_face_statistics",
    "screen_topology",
    "run_screening_batch",
    "generate_gog_configs",
    "fujii_reference_data",
    "distance_above_pareto_front",
    "plot_pareto_front",
    "plot_face_distributions",
    # general topology screening
    "compute_graph_properties",
    "screen_topology_general",
    "run_general_screening_batch",
    "generate_lattice_configs",
    "generate_random_configs",
    "plot_pc_landscape",
    "plot_pc_vs_dimension",
]


# ---------------------------------------------------------------
# Reference data
# ---------------------------------------------------------------


def fujii_reference_data() -> dict[str, dict[str, float]]:
    """Fujii & Tokunaga (2012) Table I reference thresholds.

    Returns
    -------
    dict
        Keys are geometry names, values are dicts with
        ``'p_X'``, ``'p_Z'``, ``'bpt_primal'``, ``'bpt_dual'``.
    """
    return {
        "Square": {
            "p_X": 0.103,
            "p_Z": 0.103,
            "bpt_primal": 0.500,
            "bpt_dual": 0.500,
        },
        "Kagome": {
            "p_X": 0.095,
            "p_Z": 0.116,
            "bpt_primal": 0.524,
            "bpt_dual": 0.476,
        },
        "Hexagonal": {
            "p_X": 0.065,
            "p_Z": 0.159,
            "bpt_primal": 0.653,
            "bpt_dual": 0.347,
        },
        "Tri-hexa": {
            "p_X": 0.041,
            "p_Z": 0.205,
            "bpt_primal": 0.740,
            "bpt_dual": 0.260,
        },
    }


# ---------------------------------------------------------------
# Face analysis
# ---------------------------------------------------------------


def compute_face_statistics(
    graph_factory: Callable[..., Any],
    graph_params: dict[str, Any],
    seed: int = 0,
) -> dict[str, Any]:
    """Analyse the face structure of a planar graph.

    Builds one instance, checks planarity, enumerates faces via
    :func:`~lrgsglib.graphs.nx.funcs.dual.build_planar_dual`,
    and computes face-size statistics.

    Parameters
    ----------
    graph_factory : callable
        Class returning a ``SignedGraphNX``.
    graph_params : dict
        Construction parameters (``pflip`` forced to 0).
    seed : int
        Seed for reproducibility.

    Returns
    -------
    dict
        ``'is_planar'``, ``'face_sizes'``, ``'face_size_counts'``,
        ``'uniformity_score'`` (CV of face sizes, 0 = perfectly
        uniform), ``'dominant_face_fraction'``,
        ``'n_faces'``, ``'n_vertices'``, ``'n_edges'``.
    """
    from ...graphs.nx.funcs.dual import build_planar_dual

    params = {**graph_params, "pflip": 0.0, "seed": seed}
    sg = graph_factory(**params)
    G = sg.gr[sg.on_g]

    N = G.number_of_nodes()
    E = G.number_of_edges()

    if not nx.is_planar(G):
        return {
            "is_planar": False,
            "n_vertices": N,
            "n_edges": E,
            "face_sizes": [],
            "face_size_counts": {},
            "uniformity_score": float("inf"),
            "dominant_face_fraction": 0.0,
            "n_faces": 0,
        }

    _, fi = build_planar_dual(G, include_outer_face=True)
    outer = fi["outer_face_idx"]

    sizes = [
        fi["face_sizes"][i]
        for i in range(len(fi["faces"]))
        if i != outer
    ]
    counts = Counter(sizes)
    n_faces = len(sizes)

    if n_faces == 0:
        return {
            "is_planar": True,
            "n_vertices": N,
            "n_edges": E,
            "face_sizes": sizes,
            "face_size_counts": dict(counts),
            "uniformity_score": 0.0,
            "dominant_face_fraction": 1.0,
            "n_faces": 0,
        }

    arr = np.array(sizes, dtype=float)
    mean = arr.mean()
    std = arr.std()
    cv = std / mean if mean > 0 else 0.0

    mode_size = counts.most_common(1)[0][0]
    dominant_frac = counts[mode_size] / n_faces

    return {
        "is_planar": True,
        "n_vertices": N,
        "n_edges": E,
        "face_sizes": sizes,
        "face_size_counts": dict(counts),
        "uniformity_score": float(cv),
        "dominant_face_fraction": float(dominant_frac),
        "dominant_face_size": mode_size,
        "n_faces": n_faces,
    }


# ---------------------------------------------------------------
# Topology screening
# ---------------------------------------------------------------


def screen_topology(
    graph_factory: Callable[..., Any],
    graph_params: dict[str, Any],
    p_range: NDArray[np.floating],
    n_realizations: int = 50,
    seed_base: int = 0,
    backend: str = "scipy",
    include_outer_face: bool = False,
    threshold_method: str = "fit",
    face_uniformity_cutoff: float = 0.5,
    verbose: bool = False,
) -> dict[str, Any]:
    """Screen a single topology for QEC suitability.

    Two-stage filter:

    1. **Face filter** (cheap): check planarity and face-size
       uniformity.  Non-planar or non-uniform graphs are returned
       immediately without the expensive eigenmode sweep.
    2. **Threshold computation**: run eigenmode method on primal +
       dual to extract (p_X, p_Z).

    Parameters
    ----------
    graph_factory : callable
        Class returning a ``SignedGraphNX``.
    graph_params : dict
        Construction parameters.
    p_range : ndarray
        Disorder strengths to sweep.
    n_realizations : int
        Realizations per p value.
    seed_base : int
        Base seed.
    backend : str
        Eigendecomposition backend.
    include_outer_face : bool
        Include outer face in dual.
    threshold_method : str
        ``'peak'`` or ``'fit'``.
    face_uniformity_cutoff : float
        Max coefficient of variation for faces.  Set to ``inf``
        to skip the face filter.
    verbose : bool
        Print progress.

    Returns
    -------
    dict
        ``'status'``: ``'success'``, ``'not_planar'``,
        ``'face_nonuniform'``, or ``'error'``.
        ``'face_stats'``, ``'p_X'``, ``'p_Z'``, ``'kesten_sum'``,
        ``'sweep_data'``, ``'wall_time_s'``.
    """
    t0 = time.perf_counter()

    # Stage 1: face filter
    face_stats = compute_face_statistics(
        graph_factory, graph_params, seed=seed_base
    )

    if not face_stats["is_planar"]:
        return {
            "status": "not_planar",
            "face_stats": face_stats,
            "p_X": None,
            "p_Z": None,
            "kesten_sum": None,
            "sweep_data": None,
            "wall_time_s": time.perf_counter() - t0,
        }

    if face_stats["uniformity_score"] > face_uniformity_cutoff:
        return {
            "status": "face_nonuniform",
            "face_stats": face_stats,
            "p_X": None,
            "p_Z": None,
            "kesten_sum": None,
            "sweep_data": None,
            "wall_time_s": time.perf_counter() - t0,
        }

    # Stage 2: eigenmode sweep
    try:
        p_X, p_Z, sweep = compute_pareto_point(
            graph_factory=graph_factory,
            graph_params=graph_params,
            p_range=p_range,
            n_realizations=n_realizations,
            seed_base=seed_base,
            backend=backend,
            include_outer_face=include_outer_face,
            threshold_method=threshold_method,
            verbose=verbose,
        )
    except Exception as exc:
        return {
            "status": "error",
            "face_stats": face_stats,
            "p_X": None,
            "p_Z": None,
            "kesten_sum": None,
            "sweep_data": None,
            "error": str(exc),
            "wall_time_s": time.perf_counter() - t0,
        }

    return {
        "status": "success",
        "face_stats": face_stats,
        "p_X": float(p_X),
        "p_Z": float(p_Z),
        "kesten_sum": float(p_X + p_Z),
        "sweep_data": sweep,
        "wall_time_s": time.perf_counter() - t0,
    }


# ---------------------------------------------------------------
# Batch screening
# ---------------------------------------------------------------


def run_screening_batch(
    configs: list[dict[str, Any]],
    p_range: NDArray[np.floating],
    n_realizations: int = 50,
    seed_base: int = 0,
    face_uniformity_cutoff: float = 0.5,
    verbose: bool = True,
    **sweep_kwargs: Any,
) -> pd.DataFrame:
    """Screen multiple topologies and collect results.

    Parameters
    ----------
    configs : list[dict]
        Each dict has ``'name'`` (str), ``'factory'`` (callable),
        ``'params'`` (dict).
    p_range : ndarray
        Disorder strengths.
    n_realizations : int
        Realizations per p value.
    seed_base : int
        Base seed.
    face_uniformity_cutoff : float
        Uniformity threshold for face filter.
    verbose : bool
        Print progress per config.
    **sweep_kwargs
        Forwarded to :func:`screen_topology`.

    Returns
    -------
    pd.DataFrame
        One row per config with columns: ``name``, ``status``,
        ``p_X``, ``p_Z``, ``kesten_sum``, ``uniformity_score``,
        ``dominant_face_size``, ``dominant_face_fraction``, ``N``,
        ``n_faces``, ``wall_time_s``.
    """
    rows: list[dict[str, Any]] = []

    for i, cfg in enumerate(configs):
        name = cfg["name"]
        factory = cfg["factory"]
        params = cfg["params"]

        if verbose:
            print(
                f"[{i+1}/{len(configs)}] {name}...",
                end=" ",
                flush=True,
            )

        result = screen_topology(
            graph_factory=factory,
            graph_params=params,
            p_range=p_range,
            n_realizations=n_realizations,
            seed_base=seed_base,
            face_uniformity_cutoff=face_uniformity_cutoff,
            verbose=False,
            **sweep_kwargs,
        )

        fs = result["face_stats"]
        row = {
            "name": name,
            "status": result["status"],
            "p_X": result["p_X"],
            "p_Z": result["p_Z"],
            "kesten_sum": result["kesten_sum"],
            "uniformity_score": fs.get("uniformity_score"),
            "dominant_face_size": fs.get("dominant_face_size"),
            "dominant_face_fraction": fs.get("dominant_face_fraction"),
            "N": fs.get("n_vertices"),
            "n_faces": fs.get("n_faces"),
            "wall_time_s": result["wall_time_s"],
        }
        rows.append(row)

        if verbose:
            if result["status"] == "success":
                print(
                    f"p_X={result['p_X']:.3f} p_Z={result['p_Z']:.3f} "
                    f"({result['wall_time_s']:.1f}s)"
                )
            else:
                print(f"{result['status']} ({result['wall_time_s']:.1f}s)")

    return pd.DataFrame(rows)


# ---------------------------------------------------------------
# GraphOfGraphs config generation
# ---------------------------------------------------------------


def generate_gog_configs(
    base_geos: Optional[list[str]] = None,
    fiber_geos: Optional[list[str]] = None,
    base_sides: Optional[list[int]] = None,
    fiber_sides: Optional[list[int]] = None,
    anchor_policies: Optional[list[str]] = None,
) -> list[dict[str, Any]]:
    """Generate GraphOfGraphs screening configurations.

    Enumerates the combinatorial space of backbone geometry,
    fiber geometry, sizes, and anchor policies.

    Parameters
    ----------
    base_geos : list[str]
        Base lattice geometries (default: ``['sqr', 'hex', 'tri']``).
    fiber_geos : list[str]
        Fiber lattice geometries (default: ``['sqr', 'hex', 'tri']``).
    base_sides : list[int]
        Base lattice side lengths (default: ``[3, 4, 5]``).
    fiber_sides : list[int]
        Fiber lattice side lengths (default: ``[2, 3]``).
    anchor_policies : list[str]
        Anchor policies (default: ``['first', 'center', 'random']``).

    Returns
    -------
    list[dict]
        Each dict has ``'name'``, ``'factory'``, ``'params'``.
    """
    from ...graphs.nx.GraphOfGraphsNX.GraphOfGraphsNX import (
        GraphOfGraphsNX,
    )

    if base_geos is None:
        base_geos = ["sqr", "hex", "tri"]
    if fiber_geos is None:
        fiber_geos = ["sqr", "hex", "tri"]
    if base_sides is None:
        base_sides = [3, 4, 5]
    if fiber_sides is None:
        fiber_sides = [2, 3]
    if anchor_policies is None:
        anchor_policies = ["first", "center", "random"]

    configs: list[dict[str, Any]] = []

    for bg in base_geos:
        for fg in fiber_geos:
            for bs in base_sides:
                for fs in fiber_sides:
                    for ap in anchor_policies:
                        name = f"GoG({bg}{bs}+{fg}{fs},{ap})"
                        params = {
                            "base_graph_type": "Lattice2D",
                            "base_params": {
                                "side1": bs,
                                "geo": bg,
                                "pbc": False,
                            },
                            "fiber_graph_type": "Lattice2D",
                            "fiber_params": {
                                "side1": fs,
                                "geo": fg,
                                "pbc": False,
                            },
                            "anchor_policy": ap,
                        }
                        configs.append(
                            {
                                "name": name,
                                "factory": GraphOfGraphsNX,
                                "params": params,
                            }
                        )

    return configs


# ---------------------------------------------------------------
# Pareto front analysis
# ---------------------------------------------------------------


def distance_above_pareto_front(
    p_X: float,
    p_Z: float,
    reference_points: Optional[list[tuple[float, float]]] = None,
) -> float:
    """Signed distance above the Fujii Pareto front.

    Positive = above (better than regular lattices).
    Negative = below.

    Uses perpendicular distance from the piecewise-linear
    interpolation of the reference points.

    Parameters
    ----------
    p_X, p_Z : float
        Candidate thresholds.
    reference_points : list of (p_X, p_Z) tuples, optional
        Reference Pareto front.  Defaults to Fujii data.

    Returns
    -------
    float
        Signed perpendicular distance (positive = above front).
    """
    if reference_points is None:
        fujii = fujii_reference_data()
        reference_points = sorted(
            [(d["p_X"], d["p_Z"]) for d in fujii.values()]
        )

    # Piecewise-linear interpolation of reference front
    ref = sorted(reference_points)
    ref_px = [r[0] for r in ref]
    ref_pz = [r[1] for r in ref]

    # Interpolate p_Z on the front at the candidate's p_X
    if p_X <= ref_px[0]:
        pz_front = ref_pz[0]
    elif p_X >= ref_px[-1]:
        pz_front = ref_pz[-1]
    else:
        pz_front = float(np.interp(p_X, ref_px, ref_pz))

    return p_Z - pz_front


# ---------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------


def _h_binary(p: NDArray) -> NDArray:
    """Binary Shannon entropy in bits."""
    p = np.clip(p, 1e-15, 1 - 1e-15)
    return -(p * np.log2(p) + (1 - p) * np.log2(1 - p))


def plot_pareto_front(
    results_df: pd.DataFrame,
    ax: Optional[Any] = None,
    include_gvb: bool = True,
    include_fujii: bool = True,
) -> Any:
    """Plot screening results against the Fujii Pareto front.

    Parameters
    ----------
    results_df : pd.DataFrame
        Output of :func:`run_screening_batch` (needs ``p_X``,
        ``p_Z``, ``name`` columns; only ``status='success'`` plotted).
    ax : matplotlib Axes, optional
        Target axes.  Created if None.
    include_gvb : bool
        Draw the GVB bound curve.
    include_fujii : bool
        Plot Fujii reference points.

    Returns
    -------
    matplotlib Axes
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(8, 7))

    # GVB bound
    if include_gvb:
        px = np.linspace(0.001, 0.22, 500)
        pz = []
        for pxi in px:
            target = 1.0 - _h_binary(np.array([pxi]))[0]
            if target <= 0 or target >= 1:
                pz.append(np.nan)
                continue
            lo, hi = 0.001, 0.499
            for _ in range(100):
                mid = (lo + hi) / 2
                if _h_binary(np.array([mid]))[0] < target:
                    lo = mid
                else:
                    hi = mid
            pz.append(mid)
        ax.plot(px, pz, "k--", lw=1.5, label="GVB bound")

    ax.plot(
        [0, 0.5], [0, 0.5], ":", color="gray", alpha=0.4, label="$p_X=p_Z$"
    )

    # Fujii reference
    if include_fujii:
        fujii = fujii_reference_data()
        colors = {
            "Square": "blue",
            "Kagome": "green",
            "Hexagonal": "orange",
            "Tri-hexa": "red",
        }
        for name, d in fujii.items():
            ax.plot(
                d["p_X"],
                d["p_Z"],
                "o",
                ms=10,
                color=colors.get(name, "gray"),
                label=f"{name} (Fujii)",
                zorder=5,
            )

    # Screening results
    success = results_df[results_df["status"] == "success"]
    if not success.empty:
        ax.scatter(
            success["p_X"],
            success["p_Z"],
            c="purple",
            marker="D",
            s=60,
            edgecolors="black",
            linewidths=0.5,
            zorder=6,
            label="Screened topologies",
        )
        for _, row in success.iterrows():
            ax.annotate(
                row["name"],
                (row["p_X"], row["p_Z"]),
                fontsize=6,
                textcoords="offset points",
                xytext=(5, 5),
            )

    ax.set_xlabel("$p_X$ (bit-flip threshold)", fontsize=13)
    ax.set_ylabel("$p_Z$ (phase-flip threshold)", fontsize=13)
    ax.set_title("Topology Screening: Pareto Front", fontsize=14)
    ax.legend(fontsize=8, loc="upper right")
    ax.set_xlim(0, 0.28)
    ax.set_ylim(0, 0.55)
    return ax


def plot_face_distributions(
    face_stats_list: list[dict[str, Any]],
    names: list[str],
    ax: Optional[Any] = None,
) -> Any:
    """Compare face-size distributions across topologies.

    Parameters
    ----------
    face_stats_list : list[dict]
        Each from :func:`compute_face_statistics`.
    names : list[str]
        Labels for each topology.
    ax : matplotlib Axes, optional
        Target axes.

    Returns
    -------
    matplotlib Axes
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(1, 1, figsize=(10, 5))

    all_sizes: set[int] = set()
    for fs in face_stats_list:
        all_sizes.update(fs.get("face_size_counts", {}).keys())
    sorted_sizes = sorted(all_sizes)

    n = len(names)
    width = 0.8 / n
    x = np.arange(len(sorted_sizes))

    for i, (fs, name) in enumerate(zip(face_stats_list, names)):
        counts = fs.get("face_size_counts", {})
        total = sum(counts.values()) or 1
        fracs = [counts.get(s, 0) / total for s in sorted_sizes]
        ax.bar(
            x + i * width,
            fracs,
            width,
            label=f"{name} (CV={fs.get('uniformity_score', 0):.2f})",
            alpha=0.7,
        )

    ax.set_xticks(x + width * (n - 1) / 2)
    ax.set_xticklabels(sorted_sizes)
    ax.set_xlabel("Face size (edges)")
    ax.set_ylabel("Fraction of faces")
    ax.set_title("Face-size distributions")
    ax.legend(fontsize=8)
    return ax


# ---------------------------------------------------------------
# General topology screening (any graph, not just planar)
# ---------------------------------------------------------------


def compute_graph_properties(
    graph_factory: Callable[..., Any],
    graph_params: dict[str, Any],
    seed: int = 0,
) -> dict[str, Any]:
    """Compute structural properties of a graph topology.

    Builds one instance at zero disorder and measures degree
    statistics, spectral gap, planarity, and cycle structure.

    Parameters
    ----------
    graph_factory : callable
        Class returning a ``SignedGraphNX``.
    graph_params : dict
        Construction parameters (``pflip`` forced to 0).
    seed : int
        Seed for reproducibility.

    Returns
    -------
    dict
        ``'N'``, ``'n_edges'``, ``'avg_degree'``,
        ``'degree_variance'``, ``'min_degree'``, ``'max_degree'``,
        ``'spectral_gap'``, ``'is_planar'``, ``'is_connected'``,
        ``'clustering_coefficient'``,
        ``'n_independent_cycles'`` (first Betti number b₁ = E−N+1).
    """
    params = {**graph_params, "pflip": 0.0, "seed": seed}
    sg = graph_factory(**params)
    G = sg.gr[sg.on_g]

    N = G.number_of_nodes()
    E = G.number_of_edges()
    degrees = [d for _, d in G.degree()]
    deg_arr = np.array(degrees, dtype=float)

    props: dict[str, Any] = {
        "N": N,
        "n_edges": E,
        "avg_degree": float(deg_arr.mean()) if N > 0 else 0.0,
        "degree_variance": float(deg_arr.var()) if N > 0 else 0.0,
        "min_degree": int(deg_arr.min()) if N > 0 else 0,
        "max_degree": int(deg_arr.max()) if N > 0 else 0,
        "is_planar": nx.is_planar(G),
        "is_connected": nx.is_connected(G) if N > 0 else False,
        "n_independent_cycles": max(E - N + 1, 0),
    }

    # Spectral gap (second-smallest Laplacian eigenvalue)
    if N > 2 and props["is_connected"]:
        try:
            from scipy.sparse.linalg import eigsh
            L = nx.laplacian_matrix(G).astype(float)
            vals = eigsh(L, k=2, which="SM", return_eigenvectors=False)
            props["spectral_gap"] = float(sorted(vals)[1])
        except Exception:
            eigs = np.sort(np.linalg.eigvalsh(
                nx.laplacian_matrix(G).toarray().astype(float)
            ))
            props["spectral_gap"] = float(eigs[1])
    else:
        props["spectral_gap"] = 0.0

    # Clustering coefficient
    try:
        props["clustering_coefficient"] = float(
            nx.average_clustering(G)
        )
    except Exception:
        props["clustering_coefficient"] = 0.0

    return props


def screen_topology_general(
    graph_factory: Callable[..., Any],
    graph_params: dict[str, Any],
    p_range: NDArray[np.floating],
    n_realizations: int = 50,
    seed_base: int = 0,
    backend: str = "scipy",
    threshold_method: str = "fit",
    compute_properties: bool = True,
    verbose: bool = False,
) -> dict[str, Any]:
    """Screen any topology for QEC potential via cycle dual.

    Works on **any** connected graph (planar or non-planar).
    Computes (p_X, p_Z) using the primal eigenmode and the
    minimum-cycle-basis dual eigenmode.

    Parameters
    ----------
    graph_factory : callable
        Class returning a ``SignedGraphNX``.
    graph_params : dict
        Construction parameters.
    p_range : ndarray
        Disorder strengths to sweep.
    n_realizations : int
        Realizations per p value.
    seed_base : int
        Base seed.
    backend : str
        Eigendecomposition backend.
    threshold_method : str
        ``'peak'`` or ``'fit'``.
    compute_properties : bool
        Also compute graph structural properties.
    verbose : bool
        Print progress.

    Returns
    -------
    dict
        ``'status'``, ``'p_X'``, ``'p_Z'``, ``'kesten_sum'``,
        ``'graph_properties'``, ``'sweep_data'``, ``'wall_time_s'``.
    """
    t0 = time.perf_counter()

    graph_props = None
    if compute_properties:
        graph_props = compute_graph_properties(
            graph_factory, graph_params, seed=seed_base
        )

    try:
        p_X, p_Z, sweep = compute_pareto_point_general(
            graph_factory=graph_factory,
            graph_params=graph_params,
            p_range=p_range,
            n_realizations=n_realizations,
            seed_base=seed_base,
            backend=backend,
            threshold_method=threshold_method,
            verbose=verbose,
        )
    except Exception as exc:
        return {
            "status": "error",
            "p_X": None,
            "p_Z": None,
            "kesten_sum": None,
            "graph_properties": graph_props,
            "sweep_data": None,
            "error": str(exc),
            "wall_time_s": time.perf_counter() - t0,
        }

    return {
        "status": "success",
        "p_X": float(p_X),
        "p_Z": float(p_Z),
        "kesten_sum": float(p_X + p_Z),
        "graph_properties": graph_props,
        "sweep_data": sweep,
        "wall_time_s": time.perf_counter() - t0,
    }


def run_general_screening_batch(
    configs: list[dict[str, Any]],
    p_range: NDArray[np.floating],
    n_realizations: int = 50,
    seed_base: int = 0,
    verbose: bool = True,
    **sweep_kwargs: Any,
) -> pd.DataFrame:
    """Screen multiple general topologies and collect results.

    Parameters
    ----------
    configs : list[dict]
        Each dict has ``'name'``, ``'factory'``, ``'params'``.
    p_range : ndarray
        Disorder strengths.
    n_realizations : int
        Realizations per p value.
    seed_base : int
        Base seed.
    verbose : bool
        Print progress.
    **sweep_kwargs
        Forwarded to :func:`screen_topology_general`.

    Returns
    -------
    pd.DataFrame
        Columns: ``name``, ``status``, ``p_X``, ``p_Z``,
        ``kesten_sum``, ``N``, ``avg_degree``, ``degree_variance``,
        ``spectral_gap``, ``n_cycles``, ``is_planar``,
        ``wall_time_s``.
    """
    rows: list[dict[str, Any]] = []

    for i, cfg in enumerate(configs):
        name = cfg["name"]

        if verbose:
            print(
                f"[{i+1}/{len(configs)}] {name}...",
                end=" ",
                flush=True,
            )

        result = screen_topology_general(
            graph_factory=cfg["factory"],
            graph_params=cfg["params"],
            p_range=p_range,
            n_realizations=n_realizations,
            seed_base=seed_base,
            verbose=False,
            **sweep_kwargs,
        )

        gp = result.get("graph_properties") or {}
        row = {
            "name": name,
            "status": result["status"],
            "p_X": result["p_X"],
            "p_Z": result["p_Z"],
            "kesten_sum": result["kesten_sum"],
            "N": gp.get("N"),
            "avg_degree": gp.get("avg_degree"),
            "degree_variance": gp.get("degree_variance"),
            "spectral_gap": gp.get("spectral_gap"),
            "n_cycles": gp.get("n_independent_cycles"),
            "is_planar": gp.get("is_planar"),
            "clustering": gp.get("clustering_coefficient"),
            "wall_time_s": result["wall_time_s"],
        }
        rows.append(row)

        if verbose:
            if result["status"] == "success":
                print(
                    f"p_X={result['p_X']:.3f} p_Z={result['p_Z']:.3f} "
                    f"N={gp.get('N','')} z={gp.get('avg_degree',''):.1f} "
                    f"({result['wall_time_s']:.1f}s)"
                )
            else:
                print(f"{result['status']} ({result['wall_time_s']:.1f}s)")

    return pd.DataFrame(rows)


def generate_lattice_configs(
    sizes_2d: Optional[list[int]] = None,
    sizes_3d: Optional[list[int]] = None,
    sizes_4d: Optional[list[int]] = None,
) -> list[dict[str, Any]]:
    """Generate lattice configs across dimensions (2D, 3D, 4D).

    Parameters
    ----------
    sizes_2d : list[int]
        Side lengths for 2D lattices (default: ``[16, 24]``).
    sizes_3d : list[int]
        Dimension for 3D lattices (default: ``[6, 8]``).
    sizes_4d : list[int]
        Side for 4D hypercubic (default: ``[4, 5]``).

    Returns
    -------
    list[dict]
        Each has ``'name'``, ``'factory'``, ``'params'``,
        ``'dimension'``, ``'geometry'``.
    """
    from ...graphs.nx import Lattice2DNX
    from ...graphs.nx.Lattice3DNX.Lattice3DNX import Lattice3DNX
    from ...graphs.nx.LatticeNDNX.LatticeNDNX import LatticeNDNX

    if sizes_2d is None:
        sizes_2d = [16, 24]
    if sizes_3d is None:
        sizes_3d = [6, 8]
    if sizes_4d is None:
        sizes_4d = [4, 5]

    configs: list[dict[str, Any]] = []

    # 2D lattices
    for geo, label in [("sqr", "Square"), ("hex", "Hex"), ("tri", "Tri")]:
        for s in sizes_2d:
            configs.append({
                "name": f"2D-{label}-{s}",
                "factory": Lattice2DNX,
                "params": {"side1": s, "geo": geo, "pbc": False},
                "dimension": 2,
                "geometry": geo,
            })

    # 3D lattices
    for geo, label in [("sc", "SC"), ("bcc", "BCC"), ("fcc", "FCC")]:
        for d in sizes_3d:
            configs.append({
                "name": f"3D-{label}-{d}",
                "factory": Lattice3DNX,
                "params": {"dim": d, "geo": geo, "pbc": False},
                "dimension": 3,
                "geometry": geo,
            })

    # 4D hypercubic
    for L in sizes_4d:
        configs.append({
            "name": f"4D-HC-{L}",
            "factory": LatticeNDNX,
            "params": {"shape": (L, L, L, L), "periodic": False},
            "dimension": 4,
            "geometry": "hypercubic",
        })

    return configs


def generate_random_configs(
    n_values: Optional[list[int]] = None,
    k_degrees: Optional[list[int]] = None,
) -> list[dict[str, Any]]:
    """Generate random graph configs for screening.

    Parameters
    ----------
    n_values : list[int]
        System sizes (default: ``[200, 500]``).
    k_degrees : list[int]
        Degrees for k-regular (default: ``[4, 6, 8, 10, 12, 20]``).

    Returns
    -------
    list[dict]
        Each has ``'name'``, ``'factory'``, ``'params'``,
        ``'family'``, ``'expected_avg_degree'``.
    """
    from ...graphs.nx.kRegularGraphNX.kRegularGraphNX import (
        kRegularGraphNX,
    )
    from ...graphs.nx.ErdosRenyiNX.ErdosRenyiNX import ErdosRenyiNX
    from ...graphs.nx.WattsStrogatzNX.WattsStrogatzNX import (
        WattsStrogatzNX,
    )
    from ...graphs.nx.BarabasiAlbertNX.BarabasiAlbertNX import (
        BarabasiAlbertNX,
    )

    if n_values is None:
        n_values = [200, 500]
    if k_degrees is None:
        k_degrees = [4, 6, 8, 10, 12, 20]

    configs: list[dict[str, Any]] = []

    # k-regular
    for k in k_degrees:
        for n in n_values:
            configs.append({
                "name": f"kReg-k{k}-N{n}",
                "factory": kRegularGraphNX,
                "params": {"n": n, "k": k},
                "family": "k-regular",
                "expected_avg_degree": k,
            })

    # Erdos-Renyi at matched average degrees
    for z_target in [4, 6, 8, 10]:
        for n in n_values:
            p_edge = z_target / (n - 1)
            configs.append({
                "name": f"ER-z{z_target}-N{n}",
                "factory": ErdosRenyiNX,
                "params": {"n": n, "p": p_edge},
                "family": "erdos-renyi",
                "expected_avg_degree": z_target,
            })

    # Watts-Strogatz
    for k in [4, 6, 8]:
        for p_rew in [0.01, 0.1, 0.3]:
            for n in n_values:
                configs.append({
                    "name": f"WS-k{k}-p{p_rew}-N{n}",
                    "factory": WattsStrogatzNX,
                    "params": {"n": n, "k": k, "p": p_rew},
                    "family": "watts-strogatz",
                    "expected_avg_degree": k,
                })

    # Barabasi-Albert
    for m in [2, 3, 4]:
        for n in n_values:
            configs.append({
                "name": f"BA-m{m}-N{n}",
                "factory": BarabasiAlbertNX,
                "params": {"n": n, "m": m},
                "family": "barabasi-albert",
                "expected_avg_degree": 2 * m,
            })

    return configs


# ---------------------------------------------------------------
# General topology visualization
# ---------------------------------------------------------------


def plot_pc_landscape(
    results_df: pd.DataFrame,
    x_col: str = "avg_degree",
    color_col: Optional[str] = None,
    threshold_col: str = "p_Z",
    ax: Optional[Any] = None,
) -> Any:
    """Plot threshold vs graph property, colored by family.

    Parameters
    ----------
    results_df : pd.DataFrame
        Output of :func:`run_general_screening_batch`.
    x_col : str
        Column for x-axis (e.g. ``'avg_degree'``, ``'spectral_gap'``).
    color_col : str, optional
        Column for color grouping.
    threshold_col : str
        ``'p_Z'``, ``'p_X'``, or ``'kesten_sum'``.
    ax : matplotlib Axes, optional

    Returns
    -------
    matplotlib Axes
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(9, 6))

    success = results_df[results_df["status"] == "success"].copy()
    if success.empty:
        return ax

    if color_col and color_col in success.columns:
        for label, grp in success.groupby(color_col):
            ax.scatter(
                grp[x_col], grp[threshold_col],
                label=label, s=60, edgecolors="black",
                linewidths=0.5, alpha=0.8, zorder=5,
            )
        ax.legend(fontsize=8)
    else:
        ax.scatter(
            success[x_col], success[threshold_col],
            c="steelblue", s=60, edgecolors="black",
            linewidths=0.5, zorder=5,
        )

    ax.set_xlabel(x_col.replace("_", " ").title(), fontsize=12)
    ax.set_ylabel(f"${threshold_col}$", fontsize=12)
    ax.set_title(
        f"Topology Landscape: {threshold_col} vs {x_col}", fontsize=13
    )
    return ax


def plot_pc_vs_dimension(
    results_df: pd.DataFrame,
    ax: Optional[Any] = None,
) -> Any:
    """Plot threshold vs spatial dimension for lattice families.

    Parameters
    ----------
    results_df : pd.DataFrame
        Must have ``'dimension'`` and ``'p_Z'`` columns.
    ax : matplotlib Axes, optional

    Returns
    -------
    matplotlib Axes
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 5))

    success = results_df[
        (results_df["status"] == "success")
        & results_df["dimension"].notna()
    ].copy()

    if "geometry" in success.columns:
        for geo, grp in success.groupby("geometry"):
            ax.scatter(
                grp["dimension"], grp["p_Z"],
                label=geo, s=80, edgecolors="black",
                linewidths=0.5, zorder=5,
            )
        ax.legend(fontsize=9)
    else:
        ax.scatter(
            success["dimension"], success["p_Z"],
            c="steelblue", s=80, edgecolors="black",
            linewidths=0.5, zorder=5,
        )

    ax.set_xlabel("Spatial dimension $d$", fontsize=13)
    ax.set_ylabel("$p_Z$ (Z-error threshold)", fontsize=13)
    ax.set_title("Error Threshold vs Dimension", fontsize=14)
    ax.set_xticks(sorted(success["dimension"].unique()))
    return ax
