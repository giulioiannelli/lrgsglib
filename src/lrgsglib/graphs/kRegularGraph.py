"""
Unified kRegularGraph factory with multi-engine support.

Example
-------
>>> from lrgsglib.graphs import kRegularGraph
>>>
>>> krg = kRegularGraph(n=100, k=4, engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol
    from .nx.kRegularGraphNX import kRegularGraphNX


def _get_nx_impl():
    from .nx.kRegularGraphNX import kRegularGraphNX

    return kRegularGraphNX


def _get_gt_impl():
    from .gt.kRegularGraphGT import kRegularGraphGT

    return kRegularGraphGT


register_implementation("kRegularGraph", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("kRegularGraph", GraphEngine.GRAPHTOOL, _get_gt_impl)

_NX_SPECIFIC_PARAMS = {
    "sgpathn",
    "stdFnameSFFX",
    "only_const_mode",
    "path_data",
    "path_plot",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class kRegularGraph:
    """
    Create a k-regular signed graph.

    Parameters
    ----------
    n : int
        Number of nodes.
    k : int
        Degree of each node.
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
        A k-regular graph instance.
    """

    # --- Static typing for the IDE. ---
    # This factory dispatches to a concrete engine class at runtime. We annotate
    # ``__new__`` to return the *default* engine (NetworkX) so editors give full,
    # precise method navigation -- including under ``**dict`` unpacking, which
    # defeats @overload-based typing. Params are ``Any`` to stay unpacking-proof.
    def __new__(
        cls,
        n: Any,
        k: Any,
        pflip: Any = 0.0,
        seed: Any = None,
        engine: Any = None,
        **kwargs: Any,
    ) -> "kRegularGraphNX":
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("kRegularGraph", engine)
        impl_module = impl_cls.__module__

        if "graphs.nx" in impl_module:
            impl_kwargs: dict[str, Any] = {
                "n": n,
                "k": k,
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
                "k": k,
                "pflip": pflip,
                "seed": seed,
            }
            for key in _NX_SPECIFIC_PARAMS:
                kwargs.pop(key, None)
            impl_kwargs.update(kwargs)

        else:
            impl_kwargs = {
                "n": n,
                "k": k,
                "pflip": pflip,
                "seed": seed,
                **kwargs,
            }

        return impl_cls(**impl_kwargs)
