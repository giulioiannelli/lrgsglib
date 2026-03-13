"""
Unified DualBarabasiAlbert factory with multi-engine support.

Example
-------
>>> from lrgsglib.graphs import DualBarabasiAlbert
>>>
>>> dba = DualBarabasiAlbert(n=500, m1=2, m2=1, p=0.5, engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol


def _get_nx_impl():
    from .nx.DualBarabasiAlbertNX import DualBarabasiAlbertNX

    return DualBarabasiAlbertNX


def _get_gt_impl():
    from .gt.DualBarabasiAlbertGT import DualBarabasiAlbertGT

    return DualBarabasiAlbertGT


register_implementation(
    "DualBarabasiAlbert", GraphEngine.NETWORKX, _get_nx_impl
)
register_implementation(
    "DualBarabasiAlbert", GraphEngine.GRAPHTOOL, _get_gt_impl
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


class DualBarabasiAlbert:
    """
    Create a dual Barabasi-Albert signed graph.

    Two BA processes compete: one adds m1 edges, the other m2 edges,
    chosen with probability p and 1-p respectively.

    Parameters
    ----------
    n : int
        Number of nodes.
    m1 : int
        Number of edges for the first BA process.
    m2 : int
        Number of edges for the second BA process.
    p : float, default 0.5
        Probability of using the first process.
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
        A dual BA graph instance.
    """

    def __new__(
        cls,
        n: int,
        m1: int,
        m2: int,
        p: float = 0.5,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        engine: Optional[Union[str, GraphEngine]] = None,
        **kwargs: Any,
    ):
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("DualBarabasiAlbert", engine)
        impl_module = impl_cls.__module__

        if "graphs.nx" in impl_module:
            impl_kwargs: dict[str, Any] = {
                "n": n,
                "m1": m1,
                "m2": m2,
                "p": p,
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
                "m1": m1,
                "m2": m2,
                "p": p,
                "pflip": pflip,
                "seed": seed,
            }
            for key in _NX_SPECIFIC_PARAMS:
                kwargs.pop(key, None)
            impl_kwargs.update(kwargs)

        else:
            impl_kwargs = {
                "n": n,
                "m1": m1,
                "m2": m2,
                "p": p,
                "pflip": pflip,
                "seed": seed,
                **kwargs,
            }

        return impl_cls(**impl_kwargs)
