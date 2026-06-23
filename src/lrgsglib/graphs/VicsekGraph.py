"""
Unified VicsekGraph factory with multi-engine support.

Creates Vicsek fractal graphs with hierarchical self-similar structure.

Example
-------
>>> from lrgsglib.graphs import VicsekGraph
>>>
>>> vg = VicsekGraph(N=5, k=3, engine='nx')
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

import numpy as np

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol
    from .nx.VicsekNX import VicsekGraphNX


def _get_nx_impl():
    from .nx.VicsekNX import VicsekGraphNX

    return VicsekGraphNX


def _get_gt_impl():
    from .gt.VicsekGT import VicsekGraphGT

    return VicsekGraphGT


register_implementation("VicsekGraph", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("VicsekGraph", GraphEngine.GRAPHTOOL, _get_gt_impl)

_NX_SPECIFIC_PARAMS = {
    "sgpathn",
    "stdFnameSFFX",
    "path_data",
    "path_plot",
    "init_nw_dict",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class VicsekGraph:
    """
    Create a Vicsek fractal signed graph.

    Parameters
    ----------
    N : int
        Number of motifs per level.
    k : int
        Number of hierarchical levels.
    pij : ndarray, optional
        Connection probability matrix between motifs.
    m : int, default 3
        Motif size parameter.
    symmetric : bool, default True
        Whether the motif connection matrix is symmetric.
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
        A Vicsek graph instance.
    """

    # --- Static typing for the IDE. ---
    # This factory dispatches to a concrete engine class at runtime. We annotate
    # ``__new__`` to return the *default* engine (NetworkX) so editors give full,
    # precise method navigation -- including under ``**dict`` unpacking, which
    # defeats @overload-based typing. Params are ``Any`` to stay unpacking-proof.
    def __new__(
        cls,
        N: Any,
        k: Any,
        pij: Any = None,
        m: Any = 3,
        symmetric: Any = True,
        pflip: Any = 0.0,
        seed: Any = None,
        engine: Any = None,
        **kwargs: Any,
    ) -> "VicsekGraphNX":
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("VicsekGraph", engine)
        impl_module = impl_cls.__module__

        core_kwargs: dict[str, Any] = {
            "N": N,
            "k": k,
            "m": m,
            "symmetric": symmetric,
            "pflip": pflip,
            "seed": seed,
        }
        if pij is not None:
            core_kwargs["pij"] = pij

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
