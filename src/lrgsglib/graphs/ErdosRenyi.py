"""
Unified ErdosRenyi factory with multi-engine support.

This module provides a factory function for creating Erdos-Renyi random
signed graphs using different backends (NetworkX, graph-tool, igraph).

Example
-------
>>> from lrgsglib.graphs import ErdosRenyi
>>>
>>> # Create with specific engine
>>> er_nx = ErdosRenyi(n=200, p=0.05, engine='nx')
>>> er_gt = ErdosRenyi(n=200, p=0.05, engine='gt')
>>>
>>> # Both have the same interface
>>> print(er_nx.N, er_gt.N)
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from ._base import SignedGraphProtocol


# === Lazy imports to avoid circular dependencies ===


def _get_nx_impl():
    """Lazy import for NetworkX implementation."""
    from .nx.random import ErdosRenyiNX

    return ErdosRenyiNX


def _get_gt_impl():
    """Lazy import for graph-tool implementation."""
    from .gt.random import ErdosRenyiGT

    return ErdosRenyiGT


# === Register implementations ===

register_implementation("ErdosRenyi", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("ErdosRenyi", GraphEngine.GRAPHTOOL, _get_gt_impl)


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


def ErdosRenyi(
    n: int,
    p: float,
    pflip: float = 0.0,
    extract_giant_component: bool = True,
    seed: Optional[int] = None,
    engine: Optional[Union[str, GraphEngine]] = None,
    **kwargs: Any,
) -> "SignedGraphProtocol":
    """
    Create an Erdos-Renyi random signed graph.

    This factory function creates an ER graph using the specified engine
    backend, providing a consistent interface regardless of which engine
    is used.

    Parameters
    ----------
    n : int
        Number of nodes in the initial graph.
    p : float
        Probability of edge creation (0.0 to 1.0).
    pflip : float, default 0.0
        Fraction of edges to mark for sign flipping (0.0 to 1.0).
        Call `flip_random_fract_edges()` to apply the flips.
    extract_giant_component : bool, default True
        If True, keep only the largest connected component.
        The resulting graph may have fewer nodes than n.
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
        An Erdos-Renyi graph instance supporting the signed graph protocol.

    Examples
    --------
    Create an ER graph with NetworkX:

    >>> er = ErdosRenyi(n=200, p=0.05, pflip=0.2, engine='nx')
    >>> er.flip_random_fract_edges()
    >>> print(f"Nodes: {er.N}, Negative edges: {er.count_negative_edges()}")

    Create with graph-tool (fast C++ backend):

    >>> er = ErdosRenyi(n=1000, p=0.01, engine='gt')
    >>> print(f"Nodes: {er.N}, Edges: {er.num_edges}")

    Without giant component extraction:

    >>> er = ErdosRenyi(n=100, p=0.02, extract_giant_component=False)

    Notes
    -----
    - The `pflip` parameter marks edges for flipping but doesn't apply the
      flips immediately. Call `flip_random_fract_edges()` to apply.
    - With `extract_giant_component=True`, the resulting graph will typically
      have fewer nodes than requested.
    - Both engines produce statistically equivalent graphs for the same seed.

    See Also
    --------
    lrgsglib.graphs.nx.random.ErdosRenyiNX : NetworkX implementation details
    lrgsglib.graphs.gt.random.ErdosRenyiGT : graph-tool implementation details
    """
    # Resolve engine
    if engine is not None and isinstance(engine, str):
        engine = GraphEngine(engine)

    # Get implementation class
    impl_cls = get_implementation("ErdosRenyi", engine)

    # Determine which engine we got (may have fallen back)
    impl_module = impl_cls.__module__

    # Build kwargs for the specific implementation
    if "nx_patches" in impl_module or "graphs.nx" in impl_module:
        # NetworkX implementation
        # NX ErdosRenyi always extracts GC (no option to disable)
        impl_kwargs = {
            "n": n,
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
            "p": p,
            "pflip": pflip,
            "extract_giant_component": extract_giant_component,
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
            "p": p,
            "pflip": pflip,
            "extract_giant_component": extract_giant_component,
            "seed": seed,
            **kwargs,
        }

    return impl_cls(**impl_kwargs)
