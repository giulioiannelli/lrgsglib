"""
Unified LatticeND factory with smart dimensionality dispatch.

Routes to specialized implementations when available:
- 2D shapes → Lattice2D (supports sqr/tri/hex/oct_sqr geometries)
- 3D shapes → Lattice3D (supports sc/bcc/fcc geometries)
- Any dimension → generic cubic lattice (LatticeNDNX / LatticeNDGT)

Example
-------
>>> from lrgsglib.graphs import LatticeND
>>>
>>> lat2d = LatticeND(shape=(64, 64), geo='tri')           # → Lattice2D
>>> lat3d = LatticeND(shape=(10, 10, 10), geo='bcc')       # → Lattice3D
>>> lat4d = LatticeND(shape=(5, 5, 5, 5))                  # → generic cubic
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional, Union

from ._engine import GraphEngine, get_implementation, register_implementation

if TYPE_CHECKING:
    from .protocols import SignedGraphProtocol


# === Lazy imports for generic ND ===

def _get_nx_impl():
    from .nx.LatticeNDNX import LatticeNDNX
    return LatticeNDNX


def _get_gt_impl():
    from .gt.LatticeNDGT import LatticeNDGT
    return LatticeNDGT


register_implementation("LatticeND", GraphEngine.NETWORKX, _get_nx_impl)
register_implementation("LatticeND", GraphEngine.GRAPHTOOL, _get_gt_impl)

# 2D geometries supported by Lattice2D
_GEO_2D = {"sqr", "tri", "hex", "oct_sqr"}

# 3D geometries supported by Lattice3D
_GEO_3D = {"sc", "bcc", "fcc"}

_NX_SPECIFIC_PARAMS = {
    "sgpathn",
    "stdFnameSFFX",
    "only_const_mode",
    "path_data",
    "path_plot",
}


class LatticeND:
    """
    Create an N-dimensional lattice signed graph.

    Smart dispatcher that routes to specialized implementations:

    - **2D** (``len(shape) == 2``): delegates to ``Lattice2D`` when a
      2D-specific geometry is requested (sqr, tri, hex, oct_sqr).
    - **3D** (``len(shape) == 3``): delegates to ``Lattice3D`` when a
      3D-specific geometry is requested (sc, bcc, fcc).
    - **Any dimension**: creates a generic cubic lattice using
      ``LatticeNDNX`` or ``LatticeNDGT``.

    Parameters
    ----------
    shape : tuple of int
        Dimensions of the lattice (e.g., ``(64, 64)`` for 2D,
        ``(10, 10, 10)`` for 3D, ``(5, 5, 5, 5)`` for 4D).
    geo : str, optional
        Lattice geometry. 2D: 'sqr', 'tri', 'hex', 'oct_sqr'.
        3D: 'sc', 'bcc', 'fcc'. If None, uses cubic/square default.
    periodic : bool, default True
        Whether to use periodic boundary conditions.
    pflip : float, default 0.0
        Fraction of edges to mark for sign flipping.
    seed : int, optional
        Random seed for reproducibility.
    engine : str or GraphEngine, optional
        Graph engine to use ('nx' or 'gt').
    **kwargs
        Additional engine/dimension-specific parameters (e.g.,
        ``with_positions``, ``fbc_val``, ``prew``).

    Returns
    -------
    SignedGraphProtocol
        A lattice graph instance (Lattice2D, Lattice3D, or LatticeND).

    Examples
    --------
    >>> lat = LatticeND(shape=(64, 64), geo='tri', pflip=0.3)
    >>> lat = LatticeND(shape=(10, 10, 10), geo='bcc', engine='nx')
    >>> lat = LatticeND(shape=(4, 4, 4, 4), periodic=True)
    """

    def __new__(
        cls,
        shape: tuple,
        geo: Optional[str] = None,
        periodic: bool = True,
        pflip: float = 0.0,
        seed: Optional[int] = None,
        engine: Optional[Union[str, GraphEngine]] = None,
        **kwargs: Any,
    ):
        ndim = len(shape)

        # --- Dispatch to Lattice2D for 2D shapes with 2D geometries ---
        if ndim == 2 and (geo is None or geo in _GEO_2D):
            from .Lattice2D import Lattice2D
            geo_2d = geo if geo is not None else "sqr"
            # Map generic params to Lattice2D params
            l2d_kwargs: dict[str, Any] = {
                "side1": shape[0],
                "side2": shape[1],
                "geo": geo_2d,
                "pflip": pflip,
                "seed": seed,
                "engine": engine,
            }
            # pbc/periodic mapping
            l2d_kwargs["pbc"] = periodic
            l2d_kwargs.update(kwargs)
            return Lattice2D(**l2d_kwargs)

        # --- Dispatch to Lattice3D for 3D shapes with 3D geometries ---
        if ndim == 3 and (geo is None or geo in _GEO_3D):
            from .Lattice3D import Lattice3D
            geo_3d = geo if geo is not None else "sc"
            l3d_kwargs: dict[str, Any] = {
                "dim": shape,
                "geo": geo_3d,
                "pflip": pflip,
                "seed": seed,
                "engine": engine,
            }
            # pbc/periodic mapping
            l3d_kwargs["pbc"] = periodic
            l3d_kwargs.update(kwargs)
            return Lattice3D(**l3d_kwargs)

        # --- Generic cubic lattice for any dimension ---
        if geo is not None:
            raise ValueError(
                f"Geometry '{geo}' not supported for {ndim}D lattices. "
                f"Supported: 2D={sorted(_GEO_2D)}, 3D={sorted(_GEO_3D)}. "
                f"For {ndim}D, only cubic lattice is available (geo=None)."
            )

        if engine is not None and isinstance(engine, str):
            engine = GraphEngine(engine)

        impl_cls = get_implementation("LatticeND", engine)
        impl_module = impl_cls.__module__

        impl_kwargs: dict[str, Any] = {
            "shape": shape,
            "periodic": periodic,
            "pflip": pflip,
            "seed": seed,
        }

        if "graphs.nx" in impl_module:
            for key in _NX_SPECIFIC_PARAMS:
                if key in kwargs:
                    impl_kwargs[key] = kwargs.pop(key)
            impl_kwargs.update(kwargs)
        elif "graphs.gt" in impl_module:
            for key in _NX_SPECIFIC_PARAMS:
                kwargs.pop(key, None)
            impl_kwargs.update(kwargs)
        else:
            impl_kwargs.update(kwargs)

        return impl_cls(**impl_kwargs)
