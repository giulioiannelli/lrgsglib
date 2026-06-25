"""
Unified DiracBrushGraph factory with multi-engine support.

Creates Dirac brush graphs (a graph-of-graphs with a 2D base
grid and fibers attached at each node).

Example
-------
>>> from lrgsglib.graphs import DiracBrushGraph
>>>
>>> db = DiracBrushGraph(base_x=5, base_y=5, fiber_nodes=3, engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol
    from .nx.DiracLatticeNX import DiracBrushGraphNX


def _get_nx_impl():
    from .nx.DiracLatticeNX import DiracBrushGraphNX

    return DiracBrushGraphNX


def _get_gt_impl():
    from .gt.DiracLatticeGT import DiracBrushGraphGT

    return DiracBrushGraphGT


register_implementation(
    "DiracBrushGraph", GraphEngine.NETWORKX, _get_nx_impl
)
register_implementation(
    "DiracBrushGraph", GraphEngine.GRAPHTOOL, _get_gt_impl
)

_NX_SPECIFIC_PARAMS = {
    "sgpathn",
    "stdFnameSFFX",
    "path_data",
    "path_plot",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class DiracBrushGraph:
    """
    Create a Dirac brush signed graph.

    Parameters
    ----------
    base_x : int
        X dimension of the 2D base grid.
    base_y : int
        Y dimension of the 2D base grid.
    fiber_nodes : int
        Number of nodes in each fiber.
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
        A Dirac brush graph instance.
    """

    # --- Static typing for the IDE. ---
    # Annotate ``__new__`` to return the default engine (NetworkX) so editors
    # give precise method navigation even under ``**dict`` unpacking.
    def __new__(
        cls,
        base_x: Any,
        base_y: Any,
        fiber_nodes: Any,
        periodic: Any = True,
        pflip: Any = 0.0,
        seed: Any = None,
        engine: Any = None,
        **kwargs: Any,
    ) -> "DiracBrushGraphNX":
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("DiracBrushGraph", engine)
        impl_module = impl_cls.__module__

        core_kwargs: dict[str, Any] = {
            "base_x": base_x,
            "base_y": base_y,
            "fiber_nodes": fiber_nodes,
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
