"""
Unified ExtendedBarabasiAlbert factory with multi-engine support.

Example
-------
>>> from lrgsglib.graphs import ExtendedBarabasiAlbert
>>>
>>> eba = ExtendedBarabasiAlbert(n=500, m=3, a=1.0, engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol


def _get_nx_impl():
    from .nx.ExtendedBarabasiAlbertNX import ExtendedBarabasiAlbertNX

    return ExtendedBarabasiAlbertNX


def _get_gt_impl():
    from .gt.ExtendedBarabasiAlbertGT import ExtendedBarabasiAlbertGT

    return ExtendedBarabasiAlbertGT


register_implementation(
    "ExtendedBarabasiAlbert", GraphEngine.NETWORKX, _get_nx_impl
)
register_implementation(
    "ExtendedBarabasiAlbert", GraphEngine.GRAPHTOOL, _get_gt_impl
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


class ExtendedBarabasiAlbert:
    """
    Create an extended Barabasi-Albert signed graph.

    Parameters
    ----------
    n : int
        Number of nodes.
    m : int
        Number of edges each new node creates.
    a : float, default 1.0
        Attachment exponent controlling preferential attachment strength.
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
        An extended BA graph instance.
    """

    def __new__(
        cls,
        n: int,
        m: int,
        a: float = 1.0,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        engine: Optional[Union[str, GraphEngine]] = None,
        **kwargs: Any,
    ):
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("ExtendedBarabasiAlbert", engine)
        impl_module = impl_cls.__module__

        if "graphs.nx" in impl_module:
            impl_kwargs: dict[str, Any] = {
                "n": n,
                "m": m,
                "a": a,
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
                "m": m,
                "a": a,
                "pflip": pflip,
                "seed": seed,
            }
            for key in _NX_SPECIFIC_PARAMS:
                kwargs.pop(key, None)
            impl_kwargs.update(kwargs)

        else:
            impl_kwargs = {
                "n": n,
                "m": m,
                "a": a,
                "pflip": pflip,
                "seed": seed,
                **kwargs,
            }

        return impl_cls(**impl_kwargs)
