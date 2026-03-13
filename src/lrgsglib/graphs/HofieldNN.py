"""
Unified HofieldNN factory with multi-engine support.

Example
-------
>>> from lrgsglib.graphs import HofieldNN
>>>
>>> hnn = HofieldNN(with_patterns="uniform", engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol


def _get_nx_impl():
    from .nx.HofieldNNNX import HofieldNNNX

    return HofieldNNNX


def _get_gt_impl():
    from .gt.HofieldNNGT import HofieldNNGT

    return HofieldNNGT


register_implementation("HofieldNN", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("HofieldNN", GraphEngine.GRAPHTOOL, _get_gt_impl)

_NX_SPECIFIC_PARAMS = {
    "only_const_mode",
    "path_data",
    "path_plot",
    "init_nw_dict",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class HofieldNN:
    """
    Create a Hopfield neural network signed graph.

    Parameters
    ----------
    with_patterns : str, default "uniform"
        Pattern type for Hebbian learning.
    digit : any, optional
        Specific digit/pattern to encode.
    n_samples : int, default 1
        Number of pattern samples.
    threshold : int, default 127
        Binarization threshold.
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
        A Hopfield NN graph instance.
    """

    def __new__(
        cls,
        with_patterns: str = "uniform",
        digit: Any = None,
        n_samples: int = 1,
        threshold: int = 127,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        engine: Optional[Union[str, GraphEngine]] = None,
        **kwargs: Any,
    ):
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("HofieldNN", engine)
        impl_module = impl_cls.__module__

        core_kwargs: dict[str, Any] = {
            "with_patterns": with_patterns,
            "digit": digit,
            "n_samples": n_samples,
            "threshold": threshold,
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
