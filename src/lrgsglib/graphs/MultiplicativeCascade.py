"""
Unified MultiplicativeCascade factory with multi-engine support.

Creates multiplicative cascade graphs with hierarchical structure
and tunable scale-free properties.

Example
-------
>>> from lrgsglib.graphs import MultiplicativeCascade
>>>
>>> mc = MultiplicativeCascade(p1=0.8, p2=0.6, fraction=0.4, engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol
    from .nx.MultiplicativeCascadeNX import MultiplicativeCascadeGraphNX


def _get_nx_impl():
    from .nx.MultiplicativeCascadeNX import MultiplicativeCascadeGraphNX

    return MultiplicativeCascadeGraphNX


def _get_gt_impl():
    from .gt.MultiplicativeCascadeGT import MultiplicativeCascadeGraphGT

    return MultiplicativeCascadeGraphGT


register_implementation(
    "MultiplicativeCascade", GraphEngine.NETWORKX, _get_nx_impl
)
register_implementation(
    "MultiplicativeCascade", GraphEngine.GRAPHTOOL, _get_gt_impl
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


class MultiplicativeCascade:
    """
    Create a multiplicative cascade signed graph.

    Parameters
    ----------
    p1 : float, default 0.9
        Primary cascade probability.
    p2 : float, default 0.6
        Secondary cascade probability.
    p3 : float, default 0.0
        Tertiary cascade probability.
    p4 : float, default 0.0
        Quaternary cascade probability.
    fraction : float, default 0.4
        Fraction parameter for cascade subdivision.
    iterations : int, default 5
        Number of cascade iterations.
    stochastic : bool, default False
        Use stochastic cascade variant.
    periodic : bool, default False
        Use periodic boundary conditions.
    variant : str, default "exp_clocks"
        Cascade variant type.
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
        A multiplicative cascade graph instance.
    """

    # --- Static typing for the IDE. ---
    # Annotate ``__new__`` to return the *default* engine (NetworkX) so editors
    # give full, precise method navigation -- including under ``**dict``
    # unpacking, which defeats @overload-based typing.
    def __new__(
        cls,
        p1: Any = 0.9,
        p2: Any = 0.6,
        p3: Any = 0.0,
        p4: Any = 0.0,
        fraction: Any = 0.4,
        iterations: Any = 5,
        stochastic: Any = False,
        periodic: Any = False,
        variant: Any = "exp_clocks",
        pflip: Any = 0.0,
        seed: Any = None,
        engine: Any = None,
        **kwargs: Any,
    ) -> "MultiplicativeCascadeGraphNX":
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("MultiplicativeCascade", engine)
        impl_module = impl_cls.__module__

        core_kwargs: dict[str, Any] = {
            "p1": p1,
            "p2": p2,
            "p3": p3,
            "p4": p4,
            "fraction": fraction,
            "iterations": iterations,
            "stochastic": stochastic,
            "periodic": periodic,
            "variant": variant,
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
