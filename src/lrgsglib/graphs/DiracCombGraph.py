"""
Unified DiracCombGraph factory with multi-engine support.

Creates Dirac comb graphs (a graph-of-graphs with a 1D base
and periodic fibers attached at each node).

Example
-------
>>> from lrgsglib.graphs import DiracCombGraph
>>>
>>> dc = DiracCombGraph(base_nodes=10, fiber_nodes=5, engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol


def _get_nx_impl():
    from .nx.DiracLatticeNX import DiracCombGraphNX

    return DiracCombGraphNX


def _get_gt_impl():
    from .gt.DiracLatticeGT import DiracCombGraphGT

    return DiracCombGraphGT


register_implementation(
    "DiracCombGraph", GraphEngine.NETWORKX, _get_nx_impl
)
register_implementation(
    "DiracCombGraph", GraphEngine.GRAPHTOOL, _get_gt_impl
)

_NX_SPECIFIC_PARAMS = {
    "sgpathn",
    "stdFnameSFFX",
    "path_data",
    "path_plot",
    "init_nw_dict",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class DiracCombGraph:
    """
    Create a Dirac comb signed graph.

    Parameters
    ----------
    base_nodes : int
        Number of nodes in the base graph.
    fiber_nodes : int
        Number of nodes in each fiber.
    base_type : str, default "line"
        Type of base graph.
    periodic : bool, default True
        Whether to use periodic boundary conditions.
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
        A Dirac comb graph instance.
    """

    def __new__(
        cls,
        base_nodes: int,
        fiber_nodes: int,
        base_type: str = "line",
        periodic: bool = True,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        engine: Optional[Union[str, GraphEngine]] = None,
        **kwargs: Any,
    ):
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("DiracCombGraph", engine)
        impl_module = impl_cls.__module__

        core_kwargs: dict[str, Any] = {
            "base_nodes": base_nodes,
            "fiber_nodes": fiber_nodes,
            "base_type": base_type,
            "periodic": periodic,
            "pflip": pflip,
            "seed": seed,
        }

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
