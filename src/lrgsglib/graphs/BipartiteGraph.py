"""
Unified BipartiteGraph factory with multi-engine support.

Supports both random bipartite graphs (n1, n2, p) and degree-sequence
mode (top_degrees, bottom_degrees).

Example
-------
>>> from lrgsglib.graphs import BipartiteGraph
>>>
>>> bp = BipartiteGraph(n1=50, n2=30, p=0.1, engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Sequence, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol
    from .nx.BipartiteGraphNX import BipartiteGraphNX


def _get_nx_impl():
    from .nx.BipartiteGraphNX import BipartiteGraphNX

    return BipartiteGraphNX


def _get_gt_impl():
    from .gt.BipartiteGraphGT import BipartiteGraphGT

    return BipartiteGraphGT


register_implementation(
    "BipartiteGraph", GraphEngine.NETWORKX, _get_nx_impl
)
register_implementation(
    "BipartiteGraph", GraphEngine.GRAPHTOOL, _get_gt_impl
)

_NX_SPECIFIC_PARAMS = {
    "sgpathn",
    "stdFnameSFFX",
    "only_const_mode",
    "path_data",
    "path_plot",
    "init_nw_dict",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class BipartiteGraph:
    """
    Create a bipartite signed graph.

    Supports two modes:
    - Random mode: specify n1, n2, p
    - Degree-sequence mode: specify top_degrees, bottom_degrees

    Parameters
    ----------
    n1 : int, optional
        Number of nodes in the first partition (random mode).
    n2 : int, optional
        Number of nodes in the second partition (random mode).
    p : float, default 0.1
        Edge probability between partitions (random mode).
    top_degrees : Sequence[int], optional
        Degree sequence for top partition (degree-sequence mode).
    bottom_degrees : Sequence[int], optional
        Degree sequence for bottom partition (degree-sequence mode).
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
        A bipartite graph instance.
    """

    # --- Static typing for the IDE. ---
    # This factory dispatches to a concrete engine class at runtime. We annotate
    # ``__new__`` to return the *default* engine (NetworkX) so editors give full,
    # precise method navigation -- including under ``**dict`` unpacking, which
    # defeats @overload-based typing.
    def __new__(
        cls,
        n1: Any = None,
        n2: Any = None,
        p: Any = 0.1,
        top_degrees: Any = None,
        bottom_degrees: Any = None,
        pflip: Any = 0.0,
        seed: Any = None,
        engine: Any = None,
        **kwargs: Any,
    ) -> "BipartiteGraphNX":
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("BipartiteGraph", engine)
        impl_module = impl_cls.__module__

        core_kwargs: dict[str, Any] = {
            "n1": n1,
            "n2": n2,
            "p": p,
            "pflip": pflip,
            "seed": seed,
        }
        if top_degrees is not None:
            core_kwargs["top_degrees"] = top_degrees
        if bottom_degrees is not None:
            core_kwargs["bottom_degrees"] = bottom_degrees

        if "graphs.nx" in impl_module:
            for key in _NX_SPECIFIC_PARAMS:
                if key in kwargs:
                    core_kwargs[key] = kwargs.pop(key)
            core_kwargs.update(kwargs)

        elif "graphs.gt" in impl_module:
            for key in _NX_SPECIFIC_PARAMS:
                kwargs.pop(key, None)
            core_kwargs.update(kwargs)

        else:
            core_kwargs.update(kwargs)

        return impl_cls(**core_kwargs)
