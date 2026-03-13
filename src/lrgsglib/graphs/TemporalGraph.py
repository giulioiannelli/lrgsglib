"""
Unified TemporalGraph factory with multi-engine support.

Creates temporal (time-varying) signed graphs where edges have
time intervals and signs can evolve. TemporalGraphNX inherits from
SignedGraphNX, providing full spectral/entropy/clustering functionality.

Example
-------
>>> from lrgsglib.graphs import TemporalGraph
>>>
>>> tg = TemporalGraph(n_nodes=50, pflip=0.1, seed=42, engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol


def _get_nx_impl():
    from .nx.TemporalGraphNX import TemporalGraphNX

    return TemporalGraphNX


def _get_gt_impl():
    from .gt.TemporalGraphGT import TemporalGraphGT

    return TemporalGraphGT


register_implementation(
    "TemporalGraph", GraphEngine.NETWORKX, _get_nx_impl
)
register_implementation(
    "TemporalGraph", GraphEngine.GRAPHTOOL, _get_gt_impl
)


class TemporalGraph:
    """
    Create a temporal signed graph.

    TemporalGraphNX inherits from SignedGraphNX, so all signed graph
    functionality (spectral analysis, entropy, clustering) is available.

    Parameters
    ----------
    n_nodes : int, optional
        Number of nodes.
    n_snapshots : int, optional
        Number of time snapshots.
    time_unit : str, default "continuous"
        Time model ("continuous" or "discrete").
    directed : bool, default False
        Whether edges are directed.
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
        A temporal graph instance (inherits from SignedGraphNX).
    """

    def __new__(
        cls,
        n_nodes: Optional[int] = None,
        n_snapshots: Optional[int] = None,
        time_unit: str = "continuous",
        directed: bool = False,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        engine: Optional[Union[str, GraphEngine]] = None,
        **kwargs: Any,
    ):
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("TemporalGraph", engine)

        impl_kwargs: dict[str, Any] = {
            "n_nodes": n_nodes,
            "n_snapshots": n_snapshots,
            "time_unit": time_unit,
            "directed": directed,
            "pflip": pflip,
            "seed": seed,
        }
        impl_kwargs.update(kwargs)

        return impl_cls(**impl_kwargs)
