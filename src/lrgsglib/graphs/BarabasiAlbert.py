"""
Unified BarabasiAlbert factory with multi-engine support.

This module provides a factory function for creating Barabasi-Albert
scale-free signed graphs using different backends (NetworkX, graph-tool,
igraph).

Example
-------
>>> from lrgsglib.graphs import BarabasiAlbert
>>>
>>> # Create with specific engine
>>> ba_nx = BarabasiAlbert(n=500, m=3, engine='nx')
>>> ba_gt = BarabasiAlbert(n=500, m=3, engine='gt')
>>>
>>> # Both have the same interface
>>> print(ba_nx.N, ba_gt.N)
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol
    from .nx.random import BarabasiAlbertNX


# === Lazy imports to avoid circular dependencies ===


def _get_nx_impl():
    """Lazy import for NetworkX implementation."""
    from .nx.random import BarabasiAlbertNX

    return BarabasiAlbertNX


def _get_gt_impl():
    """Lazy import for graph-tool implementation."""
    from .gt.BarabasiAlbertGT import BarabasiAlbertGT

    return BarabasiAlbertGT


# === Register implementations ===

register_implementation("BarabasiAlbert", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("BarabasiAlbert", GraphEngine.GRAPHTOOL, _get_gt_impl)


# Parameters specific to each engine (not passed to the other)
_NX_SPECIFIC_PARAMS = {
    "stdFnameSFFX",
    "sgpathn",
    "only_const_mode",
    "path_data",
    "path_plot",
    "init_nw_dict",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class BarabasiAlbert:
    """
    Barabasi-Albert scale-free signed graph.

    This class dispatches to the appropriate engine-specific implementation
    (NetworkX or graph-tool) via ``__new__``. The BA model produces
    scale-free networks with power-law degree distributions.

    Parameters
    ----------
    n : int
        Number of nodes in the final graph.
    m : int
        Number of edges each new node creates (must be <= n).
    pflip : float, default 0.0
        Fraction of edges to mark for sign flipping (0.0 to 1.0).
        Call ``flip_random_fract_edges()`` to apply the flips.
    seed : int, optional
        Random seed for reproducibility.
    engine : str or GraphEngine, optional
        Graph engine to use ('nx', 'gt', or 'ig').
        If None, uses the global default engine.
    **kwargs
        Additional engine-specific parameters.

    Examples
    --------
    >>> ba = BarabasiAlbert(n=500, m=2, pflip=0.15, engine='nx')
    >>> ba.flip_random_fract_edges()
    >>> print(f"Nodes: {ba.N}, Negative edges: {ba.count_negative_edges()}")

    See Also
    --------
    lrgsglib.graphs.nx.random.BarabasiAlbertNX : NetworkX implementation
    lrgsglib.graphs.gt.random.BarabasiAlbertGT : graph-tool implementation
    """

    # --- Static typing for the IDE. ---
    # This factory dispatches to a concrete engine class at runtime. We annotate
    # ``__new__`` to return the *default* engine (NetworkX) so editors give full,
    # precise method navigation -- including under ``**dict`` unpacking, which
    # defeats @overload-based typing.
    def __new__(
        cls,
        n: Any,
        m: Any,
        pflip: Any = 0.0,
        seed: Any = None,
        engine: Any = None,
        **kwargs: Any,
    ) -> "BarabasiAlbertNX":
        # Resolve engine
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        # Get implementation class
        impl_cls = get_implementation("BarabasiAlbert", engine)

        # Determine which engine we got (may have fallen back)
        impl_module = impl_cls.__module__

        # Build kwargs for the specific implementation
        if "graphs.nx" in impl_module:
            impl_kwargs = {
                "n": n,
                "m": m,
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
                "pflip": pflip,
                "seed": seed,
                **kwargs,
            }

        return impl_cls(**impl_kwargs)
