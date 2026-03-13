"""
Unified BipartiteFromDegreeSequence factory with multi-engine support.

Example
-------
>>> from lrgsglib.graphs import BipartiteFromDegreeSequence
>>>
>>> bp = BipartiteFromDegreeSequence(
...     top_degrees=[3, 3, 2], bottom_degrees=[2, 2, 2, 2], engine='nx'
... )
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Sequence, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol


def _get_nx_impl():
    from .nx.BipartiteFromDegreeSequenceNX import (
        BipartiteFromDegreeSequenceNX,
    )

    return BipartiteFromDegreeSequenceNX


def _get_gt_impl():
    from .gt.BipartiteFromDegreeSequenceGT import (
        BipartiteFromDegreeSequenceGT,
    )

    return BipartiteFromDegreeSequenceGT


register_implementation(
    "BipartiteFromDegreeSequence", GraphEngine.NETWORKX, _get_nx_impl
)
register_implementation(
    "BipartiteFromDegreeSequence", GraphEngine.GRAPHTOOL, _get_gt_impl
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


class BipartiteFromDegreeSequence:
    """
    Create a bipartite graph from degree sequences.

    Parameters
    ----------
    top_degrees : Sequence[int]
        Degree sequence for the top partition.
    bottom_degrees : Sequence[int]
        Degree sequence for the bottom partition.
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
        A bipartite graph instance.
    """

    def __new__(
        cls,
        top_degrees: Sequence[int],
        bottom_degrees: Sequence[int],
        pflip: float = 0.0,
        seed: Optional[int] = None,
        engine: Optional[Union[str, GraphEngine]] = None,
        **kwargs: Any,
    ):
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("BipartiteFromDegreeSequence", engine)
        impl_module = impl_cls.__module__

        core_kwargs: dict[str, Any] = {
            "top_degrees": top_degrees,
            "bottom_degrees": bottom_degrees,
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
