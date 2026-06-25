"""
Unified DGMgraph factory with multi-engine support.

Creates Dorogovtsev-Goltsev-Mendes pseudo-fractal graphs.

Example
-------
>>> from lrgsglib.graphs import DGMgraph
>>>
>>> dgm = DGMgraph(n=5, engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol
    from .nx.DGMgraphNX import DGMgraphNX


def _get_nx_impl():
    from .nx.DGMgraphNX import DGMgraphNX

    return DGMgraphNX


def _get_gt_impl():
    from .gt.DGMgraphGT import DGMgraphGT

    return DGMgraphGT


register_implementation("DGMgraph", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("DGMgraph", GraphEngine.GRAPHTOOL, _get_gt_impl)

_NX_SPECIFIC_PARAMS = {
    "sgpathn",
    "stdFnameSFFX",
    "with_positions",
    "only_const_mode",
    "path_data",
    "path_plot",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class DGMgraph:
    """
    Create a Dorogovtsev-Goltsev-Mendes pseudo-fractal signed graph.

    Parameters
    ----------
    n : int
        Fractal iteration depth.
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
        A DGM graph instance.
    """

    # --- Static typing for the IDE. ---
    # Annotate ``__new__`` to return the default engine (NetworkX) so editors
    # give precise method navigation even under ``**dict`` unpacking.
    def __new__(
        cls,
        n: Any,
        pflip: Any = 0.0,
        seed: Any = None,
        engine: Any = None,
        **kwargs: Any,
    ) -> "DGMgraphNX":
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("DGMgraph", engine)
        impl_module = impl_cls.__module__

        if "graphs.nx" in impl_module:
            impl_kwargs: dict[str, Any] = {
                "n": n,
                "pflip": pflip,
                "seed": seed,
            }
            for key in _NX_SPECIFIC_PARAMS:
                if key in kwargs:
                    impl_kwargs[key] = kwargs.pop(key)
            impl_kwargs.update(kwargs)

        elif "graphs.gt" in impl_module:
            impl_kwargs = {
                "n": n,
                "pflip": pflip,
                "seed": seed,
            }
            for key in _NX_SPECIFIC_PARAMS:
                kwargs.pop(key, None)
            impl_kwargs.update(kwargs)

        else:
            impl_kwargs = {
                "n": n,
                "pflip": pflip,
                "seed": seed,
                **kwargs,
            }

        return impl_cls(**impl_kwargs)
