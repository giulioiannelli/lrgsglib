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
    from .protocols import SignedGraphProtocol
    from .nx.random import StochasticBlockModelNX


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
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class StochasticBlockModel:
    """Create a Stochastic Block Model signed graph.

    Parameters
    ----------
    sizes : Sequence[int]
        Sizes of each community (block). Sum determines total nodes.
    p_matrix : Sequence[Sequence[float]]
        Probability matrix where p_matrix[i][j] is the edge probability
        between nodes in block i and block j. Should be symmetric.
    pflip : float, default 0.0
        Fraction of edges to mark for sign flipping (0.0 to 1.0).
    extract_giant_component : bool, default True
        If True, keep only the largest connected component.
    seed : int, optional
        Random seed for reproducibility.
    engine : str or GraphEngine, optional
        Graph engine ('nx', 'gt'). If None, uses the global default.
    **kwargs
        Additional engine-specific parameters.

    Examples
    --------
    >>> sizes = [50, 50, 50]
    >>> p_in, p_out = 0.2, 0.02
    >>> p_matrix = [[p_in, p_out, p_out],
    ...             [p_out, p_in, p_out],
    ...             [p_out, p_out, p_in]]
    >>> sbm = StochasticBlockModel(sizes=sizes, p_matrix=p_matrix, engine='nx')
    >>> sbm.flip_random_fract_edges()

    See Also
    --------
    lrgsglib.graphs.nx.random.StochasticBlockModelNX : NetworkX implementation
    lrgsglib.graphs.gt.random.StochasticBlockModelGT : graph-tool implementation
    """

    # --- Static typing for the IDE. ---
    # This factory dispatches to a concrete engine class at runtime. We annotate
    # ``__new__`` to return the *default* engine (NetworkX, the most complete
    # backend) so editors give full, precise method navigation -- including when
    # the call uses ``**dict`` unpacking, which defeats @overload-based typing.
    def __new__(
        cls,
        sizes: Any,
        p_matrix: Any,
        pflip: Any = 0.0,
        extract_giant_component: Any = True,
        seed: Any = None,
        engine: Any = None,
        **kwargs: Any,
    ) -> "StochasticBlockModelNX":
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

        elif "graphs.gt" in impl_module:
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
