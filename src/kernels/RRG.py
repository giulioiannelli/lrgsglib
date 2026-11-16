"""Random regular graph (RRG) topology for the TransCluster v2 programs.

The topology is quenched disorder too, so every rep draws its own graph:

1. Build it through the library factory with ``pflip = 0``.
2. Reduce it to canonical edges.
3. Drop the graph object.

The couplings then live in array space (``kernels/TransClusterCore.py``).
The generator seed comes from the core seed tree, so a rep always gets the
same graph whatever the tile. A disconnected draw is redrawn with the next
attempt seed.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Optional

import numpy as np
from numpy.typing import NDArray

from lrgsglib.config.progargs.defs.RandomRegular import (
    RRG_K_MIN,
    RRG_MAX_ATTEMPTS,
    RRG_PAIRING_K_MAX,
    RRG_PAIRING_MAX_EXPECTED_DRAWS,
    RRG_SAMPLER_AUTO,
    RRG_SAMPLER_NATIVE,
    RRG_SAMPLER_PAIRING,
    RRG_SAMPLERS,
)
from lrgsglib.graphs import GraphEngine, get_implementation, kRegularGraph
from lrgsglib.graphs._shared._regular_pairing import (
    expected_pairing_draws,
    validate_regular_params,
)

from .TransClusterCore import (
    canonical_edges,
    component_labels,
    edge_fingerprint,
    engine_seed,
)

__all__ = [
    "Topology",
    "draw_topology",
    "identity_engine",
    "require_engine",
    "resolve_sampler",
    "rrg_graph_key",
    "validate_rrg",
]


@dataclass(frozen=True)
class Topology:
    """Canonical edges of one rep's graph, with its provenance."""

    edges: NDArray[np.int64]
    attempt: int
    seed: int
    fingerprint: int


def validate_rrg(n_nodes: int, k: int) -> None:
    """Raise ``ValueError`` unless ``(N, k)`` is a valid RRG for this program.

    Needs ``k >= RRG_K_MIN``, since ``k = 2`` gives a union of cycles.
    """
    if k < RRG_K_MIN:
        raise ValueError(f"k must be >= {RRG_K_MIN}, got k={k}")
    validate_regular_params(n_nodes, k)


def rrg_graph_key(n_nodes: int, k: int) -> tuple[int, int]:
    """Seed-tree graph key of an RRG ensemble."""
    return (int(n_nodes), int(k))


def require_engine(engine: str) -> str:
    """``engine`` if its ``kRegularGraph`` backend imports, else
    ``ValueError``: the library would silently fall back to NX."""
    try:
        get_implementation(kRegularGraph.__name__, engine, fallback=False)
    except (ImportError, ValueError) as err:
        raise ValueError(
            f"graph engine {engine!r} is unavailable: {err}"
        ) from err
    return engine


def resolve_sampler(sampler: str, k: int, engine: str) -> str:
    """Library ``method`` for ``--sampler``.

    ``auto`` picks the exact pairing sampler up to ``RRG_PAIRING_K_MAX`` and
    the engine's native generator above it. An explicit ``pairing`` above
    that k warns, and raises once its expected draws per graph exceed
    ``RRG_PAIRING_MAX_EXPECTED_DRAWS``. Native NX is allowed, but it warns
    because it is not uniform.
    """
    if sampler not in RRG_SAMPLERS:
        raise ValueError(f"sampler must be one of {RRG_SAMPLERS}")
    method = sampler
    if sampler == RRG_SAMPLER_AUTO:
        pairing = k <= RRG_PAIRING_K_MAX
        method = RRG_SAMPLER_PAIRING if pairing else RRG_SAMPLER_NATIVE
    if method == RRG_SAMPLER_PAIRING and k > RRG_PAIRING_K_MAX:
        _check_pairing_cost(k)
    nx = GraphEngine(engine) == GraphEngine.NETWORKX
    if method == RRG_SAMPLER_NATIVE and nx:
        warnings.warn(
            "native NX random_regular_graph is not uniform (measured mean "
            "triangle count 1.61 +- 0.04 vs 4/3 at N = 1000, k = 3); do not "
            "use it for production data",
            stacklevel=2,
        )
    return method


def _check_pairing_cost(k: int) -> None:
    """Warn about the pairing cost at ``k``; raise if it is prohibitive."""
    draws = expected_pairing_draws(k)
    cost = f"the pairing sampler needs ~{draws:.3g} draws per graph at k={k}"
    if draws > RRG_PAIRING_MAX_EXPECTED_DRAWS:
        raise ValueError(f"{cost}; use --sampler native")
    warnings.warn(f"{cost}; --sampler auto would use native", stacklevel=3)


def identity_engine(method: str, engine: str) -> Optional[str]:
    """The engine, if it belongs to the ensemble identity (native only)."""
    return engine if method == RRG_SAMPLER_NATIVE else None


def draw_topology(
    n_nodes: int, k: int, *, root: int, rep: int, engine: str, method: str
) -> Topology:
    """Canonical edges of rep ``rep``'s connected k-regular graph.

    Raises ``RuntimeError`` after ``RRG_MAX_ATTEMPTS`` disconnected draws.
    """
    key = rrg_graph_key(n_nodes, k)
    for attempt in range(RRG_MAX_ATTEMPTS):
        seed = engine_seed(root, key, rep, attempt)
        graph = kRegularGraph(
            n_nodes, k, pflip=0.0, seed=seed, engine=engine, method=method
        )
        n_read, edges = canonical_edges(graph)
        del graph
        _check_regular(edges, n_read, n_nodes, k)
        if component_labels(edges, n_nodes).max() == 0:
            return Topology(edges, attempt, seed, edge_fingerprint(edges))
    raise RuntimeError(
        f"no connected {k}-regular graph on {n_nodes} nodes "
        f"in {RRG_MAX_ATTEMPTS} draws (rep {rep})"
    )


def _check_regular(
    edges: NDArray[np.int64], n_read: int, n_nodes: int, k: int
) -> None:
    """Fail loudly unless ``edges`` is a k-regular graph on ``N`` nodes."""
    degrees = np.bincount(edges.ravel(), minlength=n_nodes)
    if n_read != n_nodes or len(edges) * 2 != n_nodes * k:
        raise RuntimeError(f"generator returned {n_read} nodes, M={len(edges)}")
    if degrees.size != n_nodes or np.any(degrees != k):
        raise RuntimeError("generator returned a non-regular graph")
