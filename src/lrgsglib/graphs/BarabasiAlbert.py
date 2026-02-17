"""
Unified BarabasiAlbert factory with multi-engine support.

This module provides a factory function for creating Barabasi-Albert
scale-free signed graphs using different backends (NetworkX, graph-tool,
igraph).

Example
-------
>>> from lrgsglib.graphs import BarabasiAlbert
>>>
>>> # Create with specific engine
>>> ba_nx = BarabasiAlbert(n=500, m=3, engine='nx')
>>> ba_gt = BarabasiAlbert(n=500, m=3, engine='gt')
>>>
>>> # Both have the same interface
>>> print(ba_nx.N, ba_gt.N)
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol


# === Lazy imports to avoid circular dependencies ===


def _get_nx_impl():
    """Lazy import for NetworkX implementation."""
    from .nx.random import BarabasiAlbertNX

    return BarabasiAlbertNX


def _get_gt_impl():
    """Lazy import for graph-tool implementation."""
    from .gt.BarabasiAlbertGT import BarabasiAlbertGT

    return BarabasiAlbertGT


# === Register implementations ===

register_implementation("BarabasiAlbert", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("BarabasiAlbert", GraphEngine.GRAPHTOOL, _get_gt_impl)


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


def BarabasiAlbert(
    n: int,
    m: int,
    pflip: float = 0.0,
    seed: Optional[int] = None,
    engine: Optional[Union[str, GraphEngine]] = None,
    **kwargs: Any,
) -> "SignedGraphProtocol":
    """
    Create a Barabasi-Albert scale-free signed graph.

    This factory function creates a BA graph using the specified engine
    backend, providing a consistent interface regardless of which engine
    is used. The BA model produces scale-free networks with power-law
    degree distributions.

    Parameters
    ----------
    n : int
        Number of nodes in the final graph.
    m : int
        Number of edges each new node creates (must be <= n).
        This determines the connectivity of the resulting network.
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
        A Barabasi-Albert graph instance supporting the signed graph protocol.

    Examples
    --------
    Create a BA graph with NetworkX:

    >>> ba = BarabasiAlbert(n=500, m=2, pflip=0.15, engine='nx')
    >>> ba.flip_random_fract_edges()
    >>> print(f"Nodes: {ba.N}, Negative edges: {ba.count_negative_edges()}")

    Create with graph-tool (fast C++ backend):

    >>> ba = BarabasiAlbert(n=1000, m=3, engine='gt')
    >>> print(f"Nodes: {ba.N}, Edges: {ba.num_edges}")

    Analyze hub structure:

    >>> ba = BarabasiAlbert(n=500, m=2)
    >>> hubs = ba.get_hub_nodes(top_k=5)
    >>> print(f"Top hubs: {hubs}")

    Notes
    -----
    - The BA model produces networks with P(k) ~ k^(-gamma), gamma ~ 3.
    - The `pflip` parameter marks edges for flipping but doesn't apply the
      flips immediately. Call `flip_random_fract_edges()` to apply.
    - BA graphs are always connected by construction.
    - Expected number of edges: (n - m - 1) * m + m*(m+1)/2 for complete
      initial core.

    See Also
    --------
    lrgsglib.graphs.nx.random.BarabasiAlbertNX : NetworkX implementation
    lrgsglib.graphs.gt.random.BarabasiAlbertGT : graph-tool implementation
    """
    # Resolve engine
    if engine is not None and isinstance(engine, str):
        engine = GraphEngine(engine)

    # Get implementation class
    impl_cls = get_implementation("BarabasiAlbert", engine)

    # Determine which engine we got (may have fallen back)
    impl_module = impl_cls.__module__

    # Build kwargs for the specific implementation
    if "graphs.nx" in impl_module:
        # NetworkX implementation
        impl_kwargs = {
            "n": n,
            "m": m,
            "pflip": pflip,
            "seed": seed,
        }

        # Pass through NX-specific kwargs
        for key in _NX_SPECIFIC_PARAMS:
            if key in kwargs:
                impl_kwargs[key] = kwargs.pop(key)

        # Pass remaining kwargs
        impl_kwargs.update(kwargs)

    elif "graphs.gt" in impl_module:
        # graph-tool implementation
        impl_kwargs = {
            "n": n,
            "m": m,
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
            "m": m,
            "pflip": pflip,
            "seed": seed,
            **kwargs,
        }

    return impl_cls(**impl_kwargs)
