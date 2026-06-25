"""
Unified ErdosRenyi factory with multi-engine support.

This module provides a factory function for creating Erdos-Renyi random
signed graphs using different backends (NetworkX, graph-tool, igraph).

Example
-------
>>> from lrgsglib.graphs import ErdosRenyi
>>>
>>> # Create with specific engine
>>> er_nx = ErdosRenyi(n=200, p=0.05, engine='nx')
>>> er_gt = ErdosRenyi(n=200, p=0.05, engine='gt')
>>>
>>> # Both have the same interface
>>> print(er_nx.N, er_gt.N)
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .nx.random import ErdosRenyiNX


# === Lazy imports to avoid circular dependencies ===


def _get_nx_impl():
    """Lazy import for NetworkX implementation."""
    from .nx.random import ErdosRenyiNX

    return ErdosRenyiNX


def _get_gt_impl():
    """Lazy import for graph-tool implementation."""
    from .gt.ErdosRenyiGT import ErdosRenyiGT

    return ErdosRenyiGT


# === Register implementations ===

register_implementation("ErdosRenyi", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("ErdosRenyi", GraphEngine.GRAPHTOOL, _get_gt_impl)


# Parameters specific to each engine (not passed to the other)
_NX_SPECIFIC_PARAMS = {
    "stdFnameSFFX",
    "sgpathn",
    "only_const_mode",
    "path_data",
    "path_plot",
}

_GT_SPECIFIC_PARAMS: set[str] = set()


class ErdosRenyi:
    """Create an Erdos-Renyi random signed graph.

    Parameters
    ----------
    n : int
        Number of nodes in the initial graph.
    p : float
        Probability of edge creation (0.0 to 1.0).
    pflip : float, default 0.0
        Fraction of edges to mark for sign flipping (0.0 to 1.0).
    extract_giant_component : bool, default True
        If True, keep only the largest connected component.
    seed : int, optional
        Random seed for reproducibility.
    engine : str or GraphEngine, optional
        Graph engine ('nx', 'gt'). If None, uses the global default.
    **kwargs
        Additional engine-specific parameters.

    Examples
    --------
    >>> er = ErdosRenyi(n=200, p=0.05, pflip=0.2, engine='nx')
    >>> er.flip_random_fract_edges()

    >>> er = ErdosRenyi(n=1000, p=0.01, engine='gt')

    See Also
    --------
    lrgsglib.graphs.nx.random.ErdosRenyiNX : NetworkX implementation
    lrgsglib.graphs.gt.random.ErdosRenyiGT : graph-tool implementation
    """

    # --- Static typing for the IDE. ---
    # This factory dispatches to a concrete engine class at runtime. We annotate
    # ``__new__`` to return the *default* engine (NetworkX, the most complete
    # backend) so editors give full, precise method navigation -- including when
    # the call uses ``**dict`` unpacking, which defeats @overload-based typing
    # (pyright then falls back to the bare factory class, hiding every method).
    # For graph-tool, construct ``lrgsglib.graphs.gt.ErdosRenyiGT`` directly or
    # annotate the target.
    def __new__(
        cls,
        n: Any,
        p: Any,
        pflip: Any = 0.0,
        extract_giant_component: Any = True,
        seed: Any = None,
        engine: Any = None,
        **kwargs: Any,
    ) -> "ErdosRenyiNX":
        # Resolve engine
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        # Get implementation class
        impl_cls = get_implementation("ErdosRenyi", engine)

        # Determine which engine we got (may have fallen back)
        impl_module = impl_cls.__module__

        # Build kwargs for the specific implementation
        if "graphs.nx" in impl_module:
            # NetworkX implementation
            # NX ErdosRenyi always extracts GC (no option to disable)
            impl_kwargs = {
                "n": n,
                "p": p,
                "pflip": pflip,
                "seed": seed,
            }

            # Pass through NX-specific kwargs
            for key in _NX_SPECIFIC_PARAMS:
                if key in kwargs:
                    impl_kwargs[key] = kwargs.pop(key)

            # Pass remaining kwargs
            impl_kwargs.update(kwargs)

        elif "graphs.gt" in impl_module:
            # graph-tool implementation
            impl_kwargs = {
                "n": n,
                "p": p,
                "pflip": pflip,
                "extract_giant_component": extract_giant_component,
                "seed": seed,
            }

            # Filter out NX-specific params
            for key in _NX_SPECIFIC_PARAMS:
                kwargs.pop(key, None)

            impl_kwargs.update(kwargs)

        else:
            # Unknown implementation, pass all params
            impl_kwargs = {
                "n": n,
                "p": p,
                "pflip": pflip,
                "extract_giant_component": extract_giant_component,
                "seed": seed,
                **kwargs,
            }

        return impl_cls(**impl_kwargs)
