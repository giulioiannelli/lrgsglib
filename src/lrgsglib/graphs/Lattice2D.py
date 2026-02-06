"""
Unified Lattice2D factory with multi-engine support.

This module provides a factory function for creating 2D lattice signed graphs
using different backends (NetworkX, graph-tool, igraph).

Example
-------
>>> from lrgsglib.graphs import Lattice2D
>>>
>>> # Create with specific engine
>>> lat_nx = Lattice2D(side1=100, geo='sqr', engine='nx')
>>> lat_gt = Lattice2D(side1=100, geo='sqr', engine='gt')
>>>
>>> # Both have the same interface
>>> print(lat_nx.N, lat_gt.N)
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from ._base import SignedGraphProtocol


# === Lazy imports to avoid circular dependencies ===


def _get_nx_impl():
    """Lazy import for NetworkX implementation."""
    from .nx.Lattice2DNX import Lattice2DNX

    return Lattice2DNX


def _get_gt_impl():
    """Lazy import for graph-tool implementation."""
    from .gt.Lattice2DGT import Lattice2DGT

    return Lattice2DGT


# === Register implementations ===

register_implementation("Lattice2D", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("Lattice2D", GraphEngine.GRAPHTOOL, _get_gt_impl)


# === Parameter mapping between engines ===

# NX uses 'pbc', GT uses 'periodic'
_PARAM_MAPPING_NX_TO_GT = {
    "pbc": "periodic",
}

_PARAM_MAPPING_GT_TO_NX = {
    "periodic": "pbc",
}

# Parameters specific to each engine (not passed to the other)
_NX_SPECIFIC_PARAMS = {
    "fbc_val",
    "stdFnameSFFX",
    "sgpathn",
    "with_positions",
    "bend_positions",
    "prew",
    "only_const_mode",
}

_GT_SPECIFIC_PARAMS: set[str] = set()  # GT params are subset of common params


def Lattice2D(
    side1: int,
    side2: Optional[int] = None,
    geo: Literal["sqr", "tri", "hex"] = "sqr",
    pflip: float = 0.0,
    periodic: Optional[bool] = None,
    pbc: Optional[bool] = None,
    seed: Optional[int] = None,
    engine: Optional[Union[str, GraphEngine]] = None,
    **kwargs: Any,
) -> "SignedGraphProtocol":
    """
    Create a 2D lattice signed graph.

    This factory function creates a lattice using the specified engine backend,
    providing a consistent interface regardless of which engine is used.

    Parameters
    ----------
    side1 : int
        Size in first dimension.
    side2 : int, optional
        Size in second dimension. If None, uses side1 (square lattice).
    geo : {'sqr', 'tri', 'hex'}
        Lattice geometry:
        - 'sqr': Square lattice (4 neighbors)
        - 'tri': Triangular lattice (6 neighbors)
        - 'hex': Hexagonal lattice (3 neighbors)
    pflip : float, default 0.0
        Fraction of edges to mark for sign flipping (0.0 to 1.0).
        Call `flip_random_fract_edges()` to apply the flips.
    periodic : bool, optional
        Whether to use periodic boundary conditions.
        Alias for `pbc` (graph-tool naming convention).
    pbc : bool, optional
        Whether to use periodic boundary conditions.
        Alias for `periodic` (NetworkX naming convention).
        If both are specified, `periodic` takes precedence.
    seed : int, optional
        Random seed for reproducibility.
    engine : str or GraphEngine, optional
        Graph engine to use:
        - 'nx': NetworkX (default)
        - 'gt': graph-tool (high performance)
        - 'ig': igraph (future support)
        If None, uses the global default engine.
    **kwargs
        Additional engine-specific parameters:
        - NetworkX: fbc_val, stdFnameSFFX, sgpathn, with_positions,
          bend_positions, prew, only_const_mode, path_data, path_plot,
          init_nw_dict, backend
        - graph-tool: (currently no additional parameters)

    Returns
    -------
    SignedGraphProtocol
        A lattice graph instance supporting the signed graph protocol.

    Examples
    --------
    Create a square lattice with NetworkX:

    >>> lat = Lattice2D(side1=50, geo='sqr', pflip=0.3, engine='nx')
    >>> lat.flip_random_fract_edges()
    >>> print(f"Nodes: {lat.N}, Negative edges: {lat.count_negative_edges()}")

    Create a triangular lattice with graph-tool (fast C++ backend):

    >>> lat = Lattice2D(side1=100, geo='tri', engine='gt')
    >>> print(f"Nodes: {lat.N}, Edges: {lat.num_edges}")

    Use periodic boundary conditions:

    >>> lat = Lattice2D(side1=20, geo='sqr', periodic=True)

    Notes
    -----
    - The `pflip` parameter marks edges for flipping but doesn't apply the
      flips immediately. Call `flip_random_fract_edges()` to apply.
    - graph-tool backend offers ~50x faster graph creation for large lattices.
    - Both engines produce equivalent graph structures for the same seed.

    See Also
    --------
    lrgsglib.graphs.nx.lattice.Lattice2DNX : NetworkX implementation details
    lrgsglib.graphs.gt.lattice.Lattice2DGT : graph-tool implementation details
    """
    # Resolve engine
    if engine is not None and isinstance(engine, str):
        engine = GraphEngine(engine)

    # Get implementation class
    impl_cls = get_implementation("Lattice2D", engine)

    # Resolve periodic/pbc (periodic takes precedence)
    if periodic is None and pbc is None:
        # Default based on engine
        # NX defaults to True (L2D_PBC), GT defaults to False
        use_periodic = None  # Let implementation use its default
    elif periodic is not None:
        use_periodic = periodic
    else:
        use_periodic = pbc

    # Determine which engine we got (may have fallen back)
    impl_module = impl_cls.__module__

    # Build kwargs for the specific implementation
    if "nx_patches" in impl_module or "graphs.nx" in impl_module:
        # NetworkX implementation
        impl_kwargs = {
            "side1": side1,
            "geo": geo,
            "pflip": pflip,
            "seed": seed,
        }
        if side2 is not None:
            impl_kwargs["side2"] = side2
        if use_periodic is not None:
            impl_kwargs["pbc"] = use_periodic

        # Pass through NX-specific kwargs
        for key in _NX_SPECIFIC_PARAMS:
            if key in kwargs:
                impl_kwargs[key] = kwargs.pop(key)

        # Pass remaining kwargs (path_data, path_plot, etc.)
        impl_kwargs.update(kwargs)

    elif "gt_patches" in impl_module or "graphs.gt" in impl_module:
        # graph-tool implementation
        impl_kwargs = {
            "side1": side1,
            "side2": side2,
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
            "side1": side1,
            "side2": side2,
            "geo": geo,
            "pflip": pflip,
            "seed": seed,
            **kwargs,
        }
        if use_periodic is not None:
            impl_kwargs["periodic"] = use_periodic

    return impl_cls(**impl_kwargs)
