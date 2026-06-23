"""
Unified Lattice3D factory with multi-engine support.

This module provides a factory function for creating 3D lattice signed graphs
using different backends (NetworkX, graph-tool, igraph).

Example
-------
>>> from lrgsglib.graphs import Lattice3D
>>>
>>> # Create with specific engine
>>> lat_nx = Lattice3D(dim=10, geo='sc', engine='nx')
>>> lat_gt = Lattice3D(dim=10, geo='sc', engine='gt')
>>>
>>> # Both have the same interface
>>> print(lat_nx.N, lat_gt.N)
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .nx.lattice import Lattice3DNX


# === Lazy imports to avoid circular dependencies ===


def _get_nx_impl():
    """Lazy import for NetworkX implementation."""
    from .nx.lattice import Lattice3DNX

    return Lattice3DNX


def _get_gt_impl():
    """Lazy import for graph-tool implementation."""
    from .gt.Lattice3DGT import Lattice3DGT

    return Lattice3DGT


# === Register implementations ===

register_implementation("Lattice3D", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("Lattice3D", GraphEngine.GRAPHTOOL, _get_gt_impl)


# === Parameter mapping between engines ===

# NX uses 'pbc', GT uses 'periodic'
# NX uses 'dim' which can be int or tuple
# GT uses 'dim' similarly

# Parameters specific to NX (not passed to GT)
_NX_SPECIFIC_PARAMS = {
    "fbc_val",
    "stdFnameSFFX",
    "sgpathn",
    "with_positions",
    "pdil",
    "theta",
    "phi",
    "only_const_mode",
}


class Lattice3D:
    """Create a 3D lattice signed graph.

    Parameters
    ----------
    dim : int or tuple of 3 ints, default 10
        Lattice dimensions. If int, used for all three axes (cubic).
    geo : {'sc', 'bcc', 'fcc'}
        Lattice geometry ('sc': 6 neighbors, 'bcc': 8, 'fcc': 12).
    pflip : float, default 0.0
        Fraction of edges to mark for sign flipping (0.0 to 1.0).
    periodic / pbc : bool, optional
        Whether to use periodic boundary conditions.
        If both are specified, `periodic` takes precedence.
    seed : int, optional
        Random seed for reproducibility.
    engine : str or GraphEngine, optional
        Graph engine ('nx', 'gt'). If None, uses the global default.
    **kwargs
        Additional engine-specific parameters.

    Examples
    --------
    >>> lat = Lattice3D(dim=10, geo='sc', pflip=0.2, engine='nx')
    >>> lat.flip_random_fract_edges()

    >>> lat = Lattice3D(dim=8, geo='bcc', engine='gt', periodic=True)

    See Also
    --------
    lrgsglib.graphs.nx.lattice.Lattice3DNX : NetworkX implementation
    lrgsglib.graphs.gt.lattice.Lattice3DGT : graph-tool implementation
    """

    # --- Static typing for the IDE. ---
    # This factory dispatches to a concrete engine class at runtime. We annotate
    # ``__new__`` to return the *default* engine (NetworkX, the most complete
    # backend) so editors give full, precise method navigation -- including when
    # the call uses ``**dict`` unpacking, which defeats @overload-based typing
    # (pyright then falls back to the bare factory class, hiding every method).
    # For graph-tool, construct ``lrgsglib.graphs.gt.Lattice3DGT`` directly or
    # annotate the target.
    def __new__(
        cls,
        dim: Any = 10,
        geo: Any = "sc",
        pflip: Any = 0.0,
        periodic: Any = None,
        pbc: Any = None,
        seed: Any = None,
        engine: Any = None,
        **kwargs: Any,
    ) -> "Lattice3DNX":
        # Resolve engine
        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        # Get implementation class
        impl_cls = get_implementation("Lattice3D", engine)

        # Resolve periodic/pbc (periodic takes precedence)
        if periodic is None and pbc is None:
            use_periodic = None  # Let implementation use its default
        elif periodic is not None:
            use_periodic = periodic
        else:
            use_periodic = pbc

        # Determine which engine we got
        impl_module = impl_cls.__module__

        # Build kwargs for the specific implementation
        if "graphs.nx" in impl_module:
            # NetworkX implementation
            impl_kwargs = {
                "dim": dim,
                "geo": geo,
                "pflip": pflip,
                "seed": seed,
            }
            if use_periodic is not None:
                impl_kwargs["pbc"] = use_periodic

            # Pass through NX-specific kwargs
            for key in _NX_SPECIFIC_PARAMS:
                if key in kwargs:
                    impl_kwargs[key] = kwargs.pop(key)

            # Pass remaining kwargs
            impl_kwargs.update(kwargs)

        elif "graphs.gt" in impl_module:
            # graph-tool implementation
            impl_kwargs = {
                "dim": dim,
                "geo": geo,
                "pflip": pflip,
                "seed": seed,
            }
            if use_periodic is not None:
                impl_kwargs["periodic"] = use_periodic

            # Filter out NX-specific params
            for key in _NX_SPECIFIC_PARAMS:
                kwargs.pop(key, None)

            impl_kwargs.update(kwargs)

        else:
            # Unknown implementation, pass all params
            impl_kwargs = {
                "dim": dim,
                "geo": geo,
                "pflip": pflip,
                "seed": seed,
                **kwargs,
            }
            if use_periodic is not None:
                impl_kwargs["periodic"] = use_periodic

        return impl_cls(**impl_kwargs)
