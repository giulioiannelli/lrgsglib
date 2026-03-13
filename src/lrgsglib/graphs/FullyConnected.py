"""
Unified FullyConnected factory with multi-engine support.

Example
-------
>>> from lrgsglib.graphs import FullyConnected
>>>
>>> fc = FullyConnected(N=20, engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol


def _get_nx_impl():
    from .nx.FullyConnectedNX import FullyConnectedNX

    return FullyConnectedNX


def _get_gt_impl():
    from .gt.FullyConnectedGT import FullyConnectedGT

    return FullyConnectedGT


register_implementation(
    "FullyConnected", GraphEngine.NETWORKX, _get_nx_impl
)
register_implementation(
    "FullyConnected", GraphEngine.GRAPHTOOL, _get_gt_impl
)

_NX_SPECIFIC_PARAMS = {
    "sgpathn",
    "stdFnameSFFX",
    "only_const_mode",
    "anigemb",
    "path_data",
    "path_plot",
    "init_nw_dict",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class FullyConnected:
    """
    Create a fully connected signed graph.

    Parameters
    ----------
    N : int
        Number of nodes.
    pflip : float, default 0.0
        Fraction of edges to mark for sign flipping.
    seed : int, optional
        Random seed for reproducibility.
    with_positions : bool, default False
        Whether to compute and store node positions.
    mode_positions : str or callable, default "circular"
        Layout algorithm for node positions.
    engine : str or GraphEngine, optional
        Graph engine to use ('nx' or 'gt').
    **kwargs
        Additional engine-specific parameters.

    Returns
    -------
    SignedGraphProtocol
        A fully connected graph instance.
    """

    def __new__(
        cls,
        N: int,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        with_positions: bool = False,
        mode_positions: Union[str, Callable] = "circular",
        engine: Optional[Union[str, GraphEngine]] = None,
        **kwargs: Any,
    ):
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("FullyConnected", engine)
        impl_module = impl_cls.__module__

        if "graphs.nx" in impl_module:
            impl_kwargs: dict[str, Any] = {
                "N": N,
                "pflip": pflip,
                "seed": seed,
                "with_positions": with_positions,
                "mode_positions": mode_positions,
            }
            for key in _NX_SPECIFIC_PARAMS:
                if key in kwargs:
                    impl_kwargs[key] = kwargs.pop(key)
            impl_kwargs.update(kwargs)

        elif "graphs.gt" in impl_module:
            impl_kwargs = {
                "N": N,
                "pflip": pflip,
                "seed": seed,
                "with_positions": with_positions,
                "mode_positions": mode_positions,
            }
            for key in _NX_SPECIFIC_PARAMS:
                kwargs.pop(key, None)
            impl_kwargs.update(kwargs)

        else:
            impl_kwargs = {
                "N": N,
                "pflip": pflip,
                "seed": seed,
                "with_positions": with_positions,
                "mode_positions": mode_positions,
                **kwargs,
            }

        return impl_cls(**impl_kwargs)
