"""Helper utilities for working with :mod:`lrgsglib.Lattice3D`."""

from __future__ import annotations

from typing import Any

from lrgsglib import Lattice3D, load_or_compute_Lattice3D

__all__ = [
    "initialize_l3d_dict_args",
    "make_transcluster_lattice",
    "prepare_lattice",
]


def _normalize_dimension(dim: Any) -> Any:
    """Coerce the lattice dimension argument into the expected representation."""

    if isinstance(dim, list):
        if len(dim) == 1:
            return dim[0]
        if len(dim) == 3:
            return tuple(dim)
        raise ValueError("L must contain either 1 or 3 integers for Lattice3D")
    return dim


def initialize_l3d_dict_args(args: Any) -> dict[str, Any]:
    """Generate keyword arguments for building :class:`~lrgsglib.Lattice3D`."""

    dim = _normalize_dimension(getattr(args, "L"))
    argdict: dict[str, Any] = {
        "dim": dim,
        "geo": getattr(args, "geometry", None),
        "sgpathn": getattr(args, "workdir", None),
        "pflip": getattr(args, "p", None),
    }
    if hasattr(args, "pdil"):
        argdict["pdil"] = getattr(args, "pdil")
    return {key: value for key, value in argdict.items() if value is not None}


def make_transcluster_lattice(
    args: Any,
    *,
    init_nw_dict: bool,
    with_positions: bool = False,
    only_const_mode: bool | None = None,
) -> Lattice3D:
    """Instantiate a :class:`Lattice3D` tailored for transient cluster runs."""

    kwargs = initialize_l3d_dict_args(args)
    kwargs.update({
        "init_nw_dict": init_nw_dict,
        "with_positions": with_positions,
    })
    if only_const_mode is not None:
        kwargs["only_const_mode"] = only_const_mode
    return Lattice3D(**kwargs)


def prepare_lattice(args: Any, **kwargs: Any) -> Lattice3D:
    """Return a cached or freshly generated :class:`Lattice3D` instance."""

    base_kwargs = initialize_l3d_dict_args(args)
    cell_type = getattr(args, "cell_type", None)
    if cell_type is not None:
        base_kwargs["cell_type"] = cell_type
    return load_or_compute_Lattice3D(**base_kwargs, **kwargs)
