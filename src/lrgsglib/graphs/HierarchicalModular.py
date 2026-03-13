"""
Unified HierarchicalModular factory with multi-engine support.

Creates hierarchical modular networks with tunable community structure
at multiple scales.

Example
-------
>>> from lrgsglib.graphs import HierarchicalModular
>>>
>>> hm = HierarchicalModular(levels=3, branching=4, leaf_nodes=8, engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol


def _get_nx_impl():
    from .nx.HierarchicalModularNX import HierarchicalModularNetworkNX

    return HierarchicalModularNetworkNX


def _get_gt_impl():
    from .gt.HierarchicalModularGT import HierarchicalModularNetworkGT

    return HierarchicalModularNetworkGT


register_implementation(
    "HierarchicalModular", GraphEngine.NETWORKX, _get_nx_impl
)
register_implementation(
    "HierarchicalModular", GraphEngine.GRAPHTOOL, _get_gt_impl
)

_NX_SPECIFIC_PARAMS = {
    "sgpathn",
    "stdFnameSFFX",
    "out_suffix",
    "path_data",
    "path_plot",
    "init_nw_dict",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class HierarchicalModular:
    """
    Create a hierarchical modular signed graph.

    Parameters
    ----------
    levels : int, default 3
        Number of hierarchical levels.
    branching : int, default 4
        Branching factor at each level.
    leaf_nodes : int, default 8
        Number of nodes in each leaf module.
    p_intra : float, default 0.8
        Intra-community connection probability.
    p_ratio : float, default 0.1
        Ratio of inter/intra community edge probability per level.
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
        A hierarchical modular network instance.
    """

    def __new__(
        cls,
        levels: int = 3,
        branching: int = 4,
        leaf_nodes: int = 8,
        p_intra: float = 0.8,
        p_ratio: float = 0.1,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        engine: Optional[Union[str, GraphEngine]] = None,
        **kwargs: Any,
    ):
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("HierarchicalModular", engine)
        impl_module = impl_cls.__module__

        core_kwargs: dict[str, Any] = {
            "levels": levels,
            "branching": branching,
            "leaf_nodes": leaf_nodes,
            "p_intra": p_intra,
            "p_ratio": p_ratio,
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
