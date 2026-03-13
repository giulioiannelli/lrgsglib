"""
Unified SignedGraph factory with multi-engine support.

This module provides a factory function for creating signed graphs
using different backends (NetworkX, graph-tool).

Example
-------
>>> from lrgsglib.graphs import SignedGraph
>>>
>>> sg_nx = SignedGraph(engine='nx')
>>> sg_gt = SignedGraph(engine='gt')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol


# === Lazy imports to avoid circular dependencies ===


def _get_nx_impl():
    """Lazy import for NetworkX implementation."""
    from .nx.SignedGraphNX import SignedGraphNX

    return SignedGraphNX


def _get_gt_impl():
    """Lazy import for graph-tool implementation."""
    from .gt.SignedGraphGT import SignedGraphGT

    return SignedGraphGT


# === Register implementations ===

register_implementation("SignedGraph", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("SignedGraph", GraphEngine.GRAPHTOOL, _get_gt_impl)


# Parameters specific to each engine (not passed to the other)
_NX_SPECIFIC_PARAMS = {
    "stdFnameSFFX",
    "std_fname",
    "sgpathn",
    "only_const_mode",
    "path_data",
    "path_plot",
    "init_nw_dict",
    "init_weights_val",
    "export_mode",
    "make_dir_tree",
    "imported",
    "import_fname",
    "import_mode",
    "on_g",
    "backend",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class SignedGraph:
    """
    Create a signed graph.

    This factory class creates a signed graph using the specified engine
    backend, providing a consistent interface regardless of which engine
    is used.

    Parameters
    ----------
    G : optional
        An existing graph object to wrap as a signed graph.
        For NX: a networkx.Graph. For GT: a graph_tool.Graph.
    pflip : float, default 0.0
        Fraction of edges to mark for sign flipping.
    seed : int, optional
        Random seed for reproducibility.
    engine : str or GraphEngine, optional
        Graph engine to use ('nx' or 'gt').
    **kwargs
        Additional engine-specific parameters.

    Returns
    -------
    SignedGraphProtocol
        A signed graph instance.
    """

    def __new__(
        cls,
        G=None,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        engine: Optional[Union[str, GraphEngine]] = None,
        **kwargs: Any,
    ):
        # Resolve engine
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        # Get implementation class
        impl_cls = get_implementation("SignedGraph", engine)

        # Determine which engine we got
        impl_module = impl_cls.__module__

        # Build kwargs for the specific implementation
        if "graphs.nx" in impl_module:
            impl_kwargs: dict[str, Any] = {
                "pflip": pflip,
                "seed": seed,
            }
            if G is not None:
                impl_kwargs["G"] = G

            # Pass through NX-specific kwargs
            for key in _NX_SPECIFIC_PARAMS:
                if key in kwargs:
                    impl_kwargs[key] = kwargs.pop(key)

            impl_kwargs.update(kwargs)

        elif "graphs.gt" in impl_module:
            impl_kwargs = {
                "pflip": pflip,
                "seed": seed,
            }
            if G is not None:
                impl_kwargs["G"] = G

            # Filter out NX-specific params
            for key in _NX_SPECIFIC_PARAMS:
                kwargs.pop(key, None)

            impl_kwargs.update(kwargs)

        else:
            impl_kwargs = {
                "pflip": pflip,
                "seed": seed,
                **kwargs,
            }
            if G is not None:
                impl_kwargs["G"] = G

        return impl_cls(**impl_kwargs)
