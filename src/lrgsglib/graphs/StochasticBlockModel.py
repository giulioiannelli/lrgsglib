"""
Unified StochasticBlockModel factory with multi-engine support.

This module provides a factory function for creating Stochastic Block Model
signed graphs using different backends (NetworkX, graph-tool, igraph).

Example
-------
>>> from lrgsglib.graphs import StochasticBlockModel
>>>
>>> sizes = [30, 30, 40]
>>> p_matrix = [[0.3, 0.05, 0.05],
...             [0.05, 0.3, 0.05],
...             [0.05, 0.05, 0.3]]
>>> sbm = StochasticBlockModel(sizes=sizes, p_matrix=p_matrix, engine='gt')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Sequence, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from ._base import SignedGraphProtocol


# === Lazy imports to avoid circular dependencies ===


def _get_nx_impl():
    """Lazy import for NetworkX implementation."""
    from .nx.random import StochasticBlockModelNX

    return StochasticBlockModelNX


def _get_gt_impl():
    """Lazy import for graph-tool implementation."""
    from .gt.StochasticBlockModelGT import StochasticBlockModelGT

    return StochasticBlockModelGT


# === Register implementations ===

register_implementation("StochasticBlockModel", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("StochasticBlockModel", GraphEngine.GRAPHTOOL, _get_gt_impl)


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


def StochasticBlockModel(
    sizes: Sequence[int],
    p_matrix: Sequence[Sequence[float]],
    pflip: float = 0.0,
    extract_giant_component: bool = True,
    seed: Optional[int] = None,
    engine: Optional[Union[str, GraphEngine]] = None,
    **kwargs: Any,
) -> "SignedGraphProtocol":
    """
    Create a Stochastic Block Model signed graph.

    This factory function creates an SBM graph using the specified engine
    backend, providing a consistent interface regardless of which engine
    is used.

    Parameters
    ----------
    sizes : Sequence[int]
        Sizes of each community (block). Sum determines total nodes.
    p_matrix : Sequence[Sequence[float]]
        Probability matrix where p_matrix[i][j] is the probability
        of an edge between nodes in block i and block j.
        Should be symmetric for undirected graphs.
    pflip : float, default 0.0
        Fraction of edges to mark for sign flipping (0.0 to 1.0).
        Call `flip_random_fract_edges()` to apply the flips.
    extract_giant_component : bool, default True
        If True, keep only the largest connected component.
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
        An SBM graph instance supporting the signed graph protocol.

    Examples
    --------
    Create an SBM with three communities:

    >>> sizes = [50, 50, 50]
    >>> p_in, p_out = 0.2, 0.02
    >>> p_matrix = [[p_in, p_out, p_out],
    ...             [p_out, p_in, p_out],
    ...             [p_out, p_out, p_in]]
    >>> sbm = StochasticBlockModel(sizes=sizes, p_matrix=p_matrix, engine='nx')
    >>> sbm.flip_random_fract_edges()
    >>> print(f"Nodes: {sbm.N}, Communities: {sbm.num_communities}")

    Create with graph-tool (fast C++ backend):

    >>> sbm = StochasticBlockModel(sizes=[100, 100],
    ...                            p_matrix=[[0.3, 0.05], [0.05, 0.3]],
    ...                            engine='gt')
    >>> print(f"Modularity: {sbm.get_modularity_parameter():.3f}")

    Notes
    -----
    - The `pflip` parameter marks edges for flipping but doesn't apply the
      flips immediately. Call `flip_random_fract_edges()` to apply.
    - The SBM is useful for benchmarking community detection and studying
      spectral properties of modular networks.
    - With low inter-community probabilities, the graph may be disconnected.

    See Also
    --------
    lrgsglib.graphs.nx.random.StochasticBlockModelNX : NetworkX implementation
    lrgsglib.graphs.gt.random.StochasticBlockModelGT : graph-tool implementation
    """
    # Resolve engine
    if engine is not None and isinstance(engine, str):
        engine = GraphEngine(engine)

    # Get implementation class
    impl_cls = get_implementation("StochasticBlockModel", engine)

    # Determine which engine we got (may have fallen back)
    impl_module = impl_cls.__module__

    # Build kwargs for the specific implementation
    if "graphs.nx" in impl_module:
        # NetworkX implementation
        impl_kwargs = {
            "sizes": sizes,
            "p_matrix": p_matrix,
            "pflip": pflip,
            "seed": seed,
        }
        # NX always extracts GC (no option to disable in wrapper)

        # Pass through NX-specific kwargs
        for key in _NX_SPECIFIC_PARAMS:
            if key in kwargs:
                impl_kwargs[key] = kwargs.pop(key)

        # Pass remaining kwargs
        impl_kwargs.update(kwargs)

    elif "gt_patches" in impl_module or "graphs.gt" in impl_module:
        # graph-tool implementation
        impl_kwargs = {
            "sizes": sizes,
            "p_matrix": p_matrix,
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
            "sizes": sizes,
            "p_matrix": p_matrix,
            "pflip": pflip,
            "extract_giant_component": extract_giant_component,
            "seed": seed,
            **kwargs,
        }

    return impl_cls(**impl_kwargs)
