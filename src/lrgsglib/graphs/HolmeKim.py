"""
Unified HolmeKim factory with multi-engine support.

The Holme-Kim model extends Barabasi-Albert with triad formation
to produce scale-free networks with tunable clustering.

Example
-------
>>> from lrgsglib.graphs import HolmeKim
>>>
>>> hk = HolmeKim(n=500, m=3, p=0.5, engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol
    from .nx.HolmeKimNX import HolmeKimNX


def _get_nx_impl():
    from .nx.HolmeKimNX import HolmeKimNX

    return HolmeKimNX


def _get_gt_impl():
    from .gt.HolmeKimGT import HolmeKimGT

    return HolmeKimGT


register_implementation("HolmeKim", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("HolmeKim", GraphEngine.GRAPHTOOL, _get_gt_impl)

_NX_SPECIFIC_PARAMS = {
    "sgpathn",
    "stdFnameSFFX",
    "only_const_mode",
    "path_data",
    "path_plot",
    "init_nw_dict",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class HolmeKim:
    """
    Create a Holme-Kim signed graph (BA + clustering).

    Parameters
    ----------
    n : int
        Number of nodes.
    m : int
        Number of edges each new node creates.
    p : float, default 0.5
        Probability of triad formation step.
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
        A Holme-Kim graph instance.
    """

    # --- Static typing for the IDE. ---
    # Annotate ``__new__`` to return the *default* engine (NetworkX) so editors
    # give full, precise method navigation -- including under ``**dict``
    # unpacking, which defeats @overload-based typing.
    def __new__(
        cls,
        n: Any,
        m: Any,
        p: Any = 0.5,
        pflip: Any = 0.0,
        seed: Any = None,
        engine: Any = None,
        **kwargs: Any,
    ) -> "HolmeKimNX":
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("HolmeKim", engine)
        impl_module = impl_cls.__module__

        if "graphs.nx" in impl_module:
            impl_kwargs: dict[str, Any] = {
                "n": n,
                "m": m,
                "p": p,
                "pflip": pflip,
                "seed": seed,
            }
            for key in _NX_SPECIFIC_PARAMS:
                if key in kwargs:
                    impl_kwargs[key] = kwargs.pop(key)
            impl_kwargs.update(kwargs)

        elif "graphs.gt" in impl_module:
            impl_kwargs = {
                "n": n,
                "m": m,
                "p": p,
                "pflip": pflip,
                "seed": seed,
            }
            for key in _NX_SPECIFIC_PARAMS:
                kwargs.pop(key, None)
            impl_kwargs.update(kwargs)

        else:
            impl_kwargs = {
                "n": n,
                "m": m,
                "p": p,
                "pflip": pflip,
                "seed": seed,
                **kwargs,
            }

        return impl_cls(**impl_kwargs)
