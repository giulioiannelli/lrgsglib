"""Engine-level topology-disorder operations for the graph-tool backend.

Module-level rewiring (``prew``) and bond dilution (``pdil``) over a raw
``graph_tool.Graph``, applied by ``SignedGraphGT`` at construction time (before
the sign property is assigned) so that *every* GT graph — not just lattices —
gets these tools. The NX engine uses the parallel helpers in
``graphs/nx/funcs/base.py`` (``rewire_edges_optimized`` / ``remove_edges``); the
two implementations match each other's semantics.

Both functions are driven by ``np.random`` so a seed set by the subclass
(``np.random.seed`` + ``gt.seed_rng``) is reproducible.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from graph_tool import Graph

__all__ = ["apply_rewiring", "apply_dilution"]


def apply_rewiring(G: "Graph", prew: float) -> "Graph":
    """Rewire a fraction ``prew`` of edges in place (Watts-Strogatz style).

    One endpoint of each selected edge is reconnected to a random non-adjacent
    vertex, avoiding self-loops and duplicates. Mirrors
    ``nx.funcs.rewire_edges_optimized``. Returns the same graph object.
    """
    if prew <= 0.0:
        return G
    n_v = int(G.num_vertices())
    base_edges = [(int(u), int(v)) for u, v in G.get_edges()[:, :2]]
    existing = set()
    for u, v in base_edges:
        existing.add((u, v))
        existing.add((v, u))
    flags = np.random.rand(len(base_edges)) < prew
    for i, (u, v) in enumerate(base_edges):
        if not flags[i]:
            continue
        e = G.edge(u, v)
        if e is None:
            continue
        G.remove_edge(e)
        existing.discard((u, v))
        existing.discard((v, u))
        w = int(np.random.randint(0, n_v))
        attempts = 0
        while w == u or w == v or (u, w) in existing:
            w = int(np.random.randint(0, n_v))
            attempts += 1
            if attempts > n_v:
                break
        if w != u and G.edge(u, w) is None:
            G.add_edge(u, w)
            existing.add((u, w))
            existing.add((w, u))
    return G


def apply_dilution(
    G: "Graph", pdil: float, extract_giant: bool = True
) -> "Graph":
    """Remove a fraction ``pdil`` of edges (bond dilution).

    Dilution can disconnect the graph. When ``extract_giant`` is True the
    largest connected component is kept and returned as a fresh, compactly
    relabelled graph (vertices ``0..k-1``); vertex properties (e.g. ``pos``) are
    carried over by the filtered view. Mirrors ``nx.funcs.remove_edges``.
    """
    if pdil <= 0.0:
        return G
    import graph_tool.all as gt

    edges = G.get_edges()[:, :2]
    n_e = len(edges)
    n_remove = int(n_e * pdil)
    if n_remove <= 0:
        return G
    idx = np.random.choice(n_e, size=n_remove, replace=False)
    for u, v in edges[idx]:
        e = G.edge(int(u), int(v))
        if e is not None:
            G.remove_edge(e)

    if extract_giant:
        comp = gt.label_largest_component(G)
        if int(comp.a.sum()) < int(G.num_vertices()):
            import warnings
            warnings.warn(
                "Dilution disconnected the graph; keeping the largest "
                "connected component."
            )
            G = gt.Graph(gt.GraphView(G, vfilt=comp), prune=True)
    return G
