"""
Unified GraphOfGraphs factory with multi-engine support.

This module provides a factory function for creating hierarchical "graph of graphs"
structures using different backends (NetworkX, graph-tool).

Example
-------
>>> from lrgsglib.graphs import GraphOfGraphs
>>>
>>> # Create with specific engine
>>> gog_nx = GraphOfGraphs(
...     base_graph_type='Lattice2D', base_params={'side1': 5},
...     fiber_graph_type='ErdosRenyi', fiber_params={'n': 10, 'p': 0.2},
...     engine='nx'
... )
>>> gog_gt = GraphOfGraphs(
...     base_graph_type='Lattice2D', base_params={'side1': 5},
...     fiber_graph_type='ErdosRenyi', fiber_params={'n': 10, 'p': 0.2},
...     engine='gt'
... )
>>>
>>> # Both have the same interface
>>> print(gog_nx.N, gog_gt.N)
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol


# === Lazy imports to avoid circular dependencies ===


def _get_nx_impl():
    """Lazy import for NetworkX implementation."""
    from .nx.GraphOfGraphsNX import GraphOfGraphsNX

    return GraphOfGraphsNX


def _get_gt_impl():
    """Lazy import for graph-tool implementation."""
    from .gt.GraphOfGraphsGT import GraphOfGraphsGT

    return GraphOfGraphsGT


# === Register implementations ===

register_implementation("GraphOfGraphs", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("GraphOfGraphs", GraphEngine.GRAPHTOOL, _get_gt_impl)


class GraphOfGraphs:
    """Create a hierarchical graph-of-graphs structure.

    Each node in a base graph has a fiber graph attached to it.

    Parameters
    ----------
    base_graph_type : str
        Type name of the base graph ('Lattice2D', 'ErdosRenyi', etc.).
    base_params : dict
        Parameters to pass to the base graph constructor.
    fiber_graph_type : str
        Type name of the fiber graph.
    fiber_params : dict or callable
        Parameters for fiber graph construction.
        If callable: function(base_idx) -> dict for heterogeneous fibers.
    anchor_policy : str or callable, default 'first'
        How to select anchor nodes ('first', 'center', 'last', 'random',
        or callable(base_idx, n_fiber) -> anchor_idx).
    pflip : float, default 0.0
        Fraction of edges to mark for sign flipping (0.0 to 1.0).
    seed : int, optional
        Random seed for reproducibility.
    engine : str or GraphEngine, optional
        Graph engine ('nx', 'gt'). If None, uses the global default.
    **kwargs
        Additional engine-specific parameters.

    Examples
    --------
    >>> gog = GraphOfGraphs(
    ...     base_graph_type='Lattice2D',
    ...     base_params={'side1': 5, 'geo': 'sqr'},
    ...     fiber_graph_type='ErdosRenyi',
    ...     fiber_params={'n': 10, 'p': 0.2},
    ...     engine='nx'
    ... )

    See Also
    --------
    lrgsglib.graphs.nx.GraphOfGraphsNX : NetworkX implementation
    lrgsglib.graphs.gt.GraphOfGraphsGT : graph-tool implementation
    """

    def __new__(
        cls,
        base_graph_type: str,
        base_params: dict,
        fiber_graph_type: str,
        fiber_params: Union[dict, Callable[[int], dict]],
        anchor_policy: Union[str, Callable[[int, int], int]] = "first",
        pflip: float = 0.0,
        seed: Optional[int] = None,
        engine: Optional[Union[str, GraphEngine]] = None,
        **kwargs: Any,
    ):
        # Resolve engine
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        # Get implementation class
        impl_cls = get_implementation("GraphOfGraphs", engine)

        # Build kwargs for the implementation
        impl_kwargs = {
            "base_graph_type": base_graph_type,
            "base_params": base_params,
            "fiber_graph_type": fiber_graph_type,
            "fiber_params": fiber_params,
            "anchor_policy": anchor_policy,
            "pflip": pflip,
            "seed": seed,
        }

        # Pass through additional kwargs
        impl_kwargs.update(kwargs)

        return impl_cls(**impl_kwargs)
