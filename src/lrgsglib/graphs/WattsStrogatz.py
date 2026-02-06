"""
Unified WattsStrogatz factory with multi-engine support.

This module provides a factory function for creating Watts-Strogatz small-world
signed graphs using different backends (NetworkX, graph-tool, igraph).

Example
-------
>>> from lrgsglib.graphs import WattsStrogatz
>>>
>>> # Create with specific engine
>>> ws_nx = WattsStrogatz(n=100, k=4, p=0.3, engine='nx')
>>> ws_gt = WattsStrogatz(n=100, k=4, p=0.3, engine='gt')
>>>
>>> # Both have the same interface
>>> print(ws_nx.N, ws_gt.N)
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from ._base import SignedGraphProtocol


# === Lazy imports to avoid circular dependencies ===


def _get_nx_impl():
    """Lazy import for NetworkX implementation."""
    from .nx.random import WattsStrogatzNX

    return WattsStrogatzNX


def _get_gt_impl():
    """Lazy import for graph-tool implementation."""
    from .gt.WattsStrogatzGT import WattsStrogatzGT

    return WattsStrogatzGT


# === Register implementations ===

register_implementation("WattsStrogatz", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("WattsStrogatz", GraphEngine.GRAPHTOOL, _get_gt_impl)


# Parameters specific to each engine (not passed to the other)
_NX_SPECIFIC_PARAMS = {
    "stdFnameSFFX",
    "sgpathn",
    "only_const_mode",
    "path_data",
    "path_plot",
    "init_nw_dict",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


def WattsStrogatz(
    n: int,
    k: int,
    p: float,
    pflip: float = 0.0,
    seed: Optional[int] = None,
    engine: Optional[Union[str, GraphEngine]] = None,
    **kwargs: Any,
) -> "SignedGraphProtocol":
    """
    Create a Watts-Strogatz small-world signed graph.

    This factory function creates a WS graph using the specified engine
    backend, providing a consistent interface regardless of which engine
    is used.

    Parameters
    ----------
    n : int
        Number of nodes in the graph.
    k : int
        Each node is connected to k nearest neighbors in ring topology.
        Must be even.
    p : float
        Probability of rewiring each edge (0.0 to 1.0).
    pflip : float, default 0.0
        Fraction of edges to mark for sign flipping (0.0 to 1.0).
        Call `flip_random_fract_edges()` to apply the flips.
    seed : int, optional
        Random seed for reproducibility.
    engine : str or GraphEngine, optional
        Graph engine to use:
        - 'nx': NetworkX (default)
        - 'gt': graph-tool (high performance)
        - 'ig': igraph (future support)
        If None, uses the global default engine.
    **kwargs
        Additional engine-specific parameters:
        - NetworkX: stdFnameSFFX, sgpathn, only_const_mode, path_data,
          path_plot, init_nw_dict
        - graph-tool: (currently no additional parameters)

    Returns
    -------
    SignedGraphProtocol
        A Watts-Strogatz graph instance supporting the signed graph protocol.

    Examples
    --------
    Create a WS graph with NetworkX:

    >>> ws = WattsStrogatz(n=100, k=4, p=0.3, pflip=0.1, engine='nx')
    >>> ws.flip_random_fract_edges()
    >>> print(f"Nodes: {ws.N}, Negative edges: {ws.count_negative_edges()}")

    Create with graph-tool (fast C++ backend):

    >>> ws = WattsStrogatz(n=1000, k=6, p=0.2, engine='gt')
    >>> print(f"Nodes: {ws.N}, Edges: {ws.num_edges}")

    Notes
    -----
    - The `pflip` parameter marks edges for flipping but doesn't apply the
      flips immediately. Call `flip_random_fract_edges()` to apply.
    - The WS model interpolates between regular lattices (p=0) and random
      graphs (p=1), producing small-world networks for intermediate p.
    - NetworkX uses connected_watts_strogatz_graph (always connected).

    See Also
    --------
    lrgsglib.graphs.nx.random.WattsStrogatzNX : NetworkX implementation
    lrgsglib.graphs.gt.random.WattsStrogatzGT : graph-tool implementation
    """
    # Resolve engine
    if engine is not None and isinstance(engine, str):
        engine = GraphEngine(engine)

    # Get implementation class
    impl_cls = get_implementation("WattsStrogatz", engine)

    # Determine which engine we got (may have fallen back)
    impl_module = impl_cls.__module__

    # Build kwargs for the specific implementation
    if "nx_patches" in impl_module or "graphs.nx" in impl_module:
        # NetworkX implementation
        impl_kwargs = {
            "n": n,
            "k": k,
            "p": p,
            "pflip": pflip,
            "seed": seed,
        }

        # Pass through NX-specific kwargs
        for key in _NX_SPECIFIC_PARAMS:
            if key in kwargs:
                impl_kwargs[key] = kwargs.pop(key)

        # Pass remaining kwargs
        impl_kwargs.update(kwargs)

    elif "gt_patches" in impl_module or "graphs.gt" in impl_module:
        # graph-tool implementation
        impl_kwargs = {
            "n": n,
            "k": k,
            "p": p,
            "pflip": pflip,
            "seed": seed,
        }

        # Filter out NX-specific params
        for key in _NX_SPECIFIC_PARAMS:
            kwargs.pop(key, None)

        impl_kwargs.update(kwargs)

    else:
        # Unknown implementation, pass all params
        impl_kwargs = {
            "n": n,
            "k": k,
            "p": p,
            "pflip": pflip,
            "seed": seed,
            **kwargs,
        }

    return impl_cls(**impl_kwargs)
