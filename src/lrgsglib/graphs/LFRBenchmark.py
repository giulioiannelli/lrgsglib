"""
Unified LFRBenchmark factory with multi-engine support.

The Lancichinetti-Fortunato-Radicchi benchmark generates networks
with planted community structure and power-law degree distributions.

Example
-------
>>> from lrgsglib.graphs import LFRBenchmark
>>>
>>> lfr = LFRBenchmark(n=250, tau1=3, tau2=1.5, mu=0.1, engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol


def _get_nx_impl():
    from .nx.LFRBenchmarkNX import LFRBenchmarkNX

    return LFRBenchmarkNX


def _get_gt_impl():
    from .gt.LFRBenchmarkGT import LFRBenchmarkGT

    return LFRBenchmarkGT


register_implementation("LFRBenchmark", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("LFRBenchmark", GraphEngine.GRAPHTOOL, _get_gt_impl)

_NX_SPECIFIC_PARAMS = {
    "sgpathn",
    "stdFnameSFFX",
    "only_const_mode",
    "path_data",
    "path_plot",
    "init_nw_dict",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class LFRBenchmark:
    """
    Create an LFR benchmark signed graph.

    Parameters
    ----------
    n : int
        Number of nodes.
    tau1 : float
        Power-law exponent for degree distribution.
    tau2 : float
        Power-law exponent for community size distribution.
    mu : float
        Mixing parameter (fraction of inter-community edges).
    average_degree : float, optional
        Average node degree.
    min_degree : int, optional
        Minimum node degree.
    max_degree : int, optional
        Maximum node degree.
    min_community : int, optional
        Minimum community size.
    max_community : int, optional
        Maximum community size.
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
        An LFR benchmark graph instance.
    """

    def __new__(
        cls,
        n: int,
        tau1: float,
        tau2: float,
        mu: float,
        average_degree: Optional[float] = None,
        min_degree: Optional[int] = None,
        max_degree: Optional[int] = None,
        min_community: Optional[int] = None,
        max_community: Optional[int] = None,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        engine: Optional[Union[str, GraphEngine]] = None,
        **kwargs: Any,
    ):
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("LFRBenchmark", engine)
        impl_module = impl_cls.__module__

        core_kwargs: dict[str, Any] = {
            "n": n,
            "tau1": tau1,
            "tau2": tau2,
            "mu": mu,
            "pflip": pflip,
            "seed": seed,
        }
        if average_degree is not None:
            core_kwargs["average_degree"] = average_degree
        if min_degree is not None:
            core_kwargs["min_degree"] = min_degree
        if max_degree is not None:
            core_kwargs["max_degree"] = max_degree
        if min_community is not None:
            core_kwargs["min_community"] = min_community
        if max_community is not None:
            core_kwargs["max_community"] = max_community

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
