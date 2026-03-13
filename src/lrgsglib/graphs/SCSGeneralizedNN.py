"""
Unified SCSGeneralizedNN factory with multi-engine support.

Example
-------
>>> from lrgsglib.graphs import SCSGeneralizedNN
>>>
>>> scs = SCSGeneralizedNN(N=100, gamma=0.5, engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol


def _get_nx_impl():
    from .nx.SCSGeneralizedNNNX import SCSGeneralizedNNNX

    return SCSGeneralizedNNNX


def _get_gt_impl():
    from .gt.SCSGeneralizedNNGT import SCSGeneralizedNNGT

    return SCSGeneralizedNNGT


register_implementation(
    "SCSGeneralizedNN", GraphEngine.NETWORKX, _get_nx_impl
)
register_implementation(
    "SCSGeneralizedNN", GraphEngine.GRAPHTOOL, _get_gt_impl
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


class SCSGeneralizedNN:
    """
    Create an SCS generalized neural network signed graph.

    Parameters
    ----------
    N : int
        Number of nodes.
    gamma : float
        Coupling strength parameter.
    J0 : float, default 0.0
        Mean coupling.
    J : float, default 1.0
        Coupling variance.
    g : float, default 1.0
        Gain parameter.
    diagonal : any, optional
        Diagonal elements of the coupling matrix.
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
        An SCS generalized NN graph instance.
    """

    def __new__(
        cls,
        N: int,
        gamma: float,
        J0: float = 0.0,
        J: float = 1.0,
        g: float = 1.0,
        diagonal: Any = None,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        engine: Optional[Union[str, GraphEngine]] = None,
        **kwargs: Any,
    ):
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("SCSGeneralizedNN", engine)
        impl_module = impl_cls.__module__

        core_kwargs: dict[str, Any] = {
            "N": N,
            "gamma": gamma,
            "J0": J0,
            "J": J,
            "g": g,
            "pflip": pflip,
            "seed": seed,
        }
        if diagonal is not None:
            core_kwargs["diagonal"] = diagonal

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
