"""
Unified ConfigurationModel factory with multi-engine support.

Example
-------
>>> from lrgsglib.graphs import ConfigurationModel
>>>
>>> cm = ConfigurationModel(degree_sequence=[3, 3, 3, 3, 2, 2], engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Sequence, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol
    from .nx.ConfigurationModelNX import ConfigurationModelNX


def _get_nx_impl():
    from .nx.ConfigurationModelNX import ConfigurationModelNX

    return ConfigurationModelNX


def _get_gt_impl():
    from .gt.ConfigurationModelGT import ConfigurationModelGT

    return ConfigurationModelGT


register_implementation(
    "ConfigurationModel", GraphEngine.NETWORKX, _get_nx_impl
)
register_implementation(
    "ConfigurationModel", GraphEngine.GRAPHTOOL, _get_gt_impl
)

_NX_SPECIFIC_PARAMS = {
    "sgpathn",
    "stdFnameSFFX",
    "only_const_mode",
    "path_data",
    "path_plot",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class ConfigurationModel:
    """
    Create a graph from a degree sequence using the configuration model.

    Parameters
    ----------
    degree_sequence : Sequence[int]
        The degree sequence to realize.
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
        A configuration model graph instance.
    """

    # --- Static typing for the IDE. ---
    # This factory dispatches to a concrete engine class at runtime. We annotate
    # ``__new__`` to return the *default* engine (NetworkX) so editors give full,
    # precise method navigation -- including under ``**dict`` unpacking, which
    # defeats @overload-based typing.
    def __new__(
        cls,
        degree_sequence: Any,
        pflip: Any = 0.0,
        seed: Any = None,
        engine: Any = None,
        **kwargs: Any,
    ) -> "ConfigurationModelNX":
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("ConfigurationModel", engine)
        impl_module = impl_cls.__module__

        if "graphs.nx" in impl_module:
            impl_kwargs: dict[str, Any] = {
                "degree_sequence": degree_sequence,
                "pflip": pflip,
                "seed": seed,
            }
            for key in _NX_SPECIFIC_PARAMS:
                if key in kwargs:
                    impl_kwargs[key] = kwargs.pop(key)
            impl_kwargs.update(kwargs)

        elif "graphs.gt" in impl_module:
            impl_kwargs = {
                "degree_sequence": degree_sequence,
                "pflip": pflip,
                "seed": seed,
            }
            for key in _NX_SPECIFIC_PARAMS:
                kwargs.pop(key, None)
            impl_kwargs.update(kwargs)

        else:
            impl_kwargs = {
                "degree_sequence": degree_sequence,
                "pflip": pflip,
                "seed": seed,
                **kwargs,
            }

        return impl_cls(**impl_kwargs)
