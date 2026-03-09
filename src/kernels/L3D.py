"""Helper utilities for working with :mod:`lrgsglib.Lattice3D`."""

from __future__ import annotations

from typing import Any

from lrgsglib import Lattice3D, load_or_compute_Lattice3D
from lrgsglib.utils.basic.strings import join_non_empty
from parsers.shared import resolve_graph_class, get_graph_engine

__all__ = [
    "get_l3d_out_suffix",
    "initialize_l3d_dict_args",
    "make_transcluster_lattice",
    "prepare_lattice",
]


def get_l3d_out_suffix(args: Any) -> str:
    """Build output suffix from L3D-relevant structural parameters."""

    init_cond = getattr(args, "init_cond", "")
    geometry = getattr(args, "geometry", "")
    pdil = getattr(args, "pdil", 0.0)
    edge_weight = getattr(args, "edge_weight", "flip")
    mu = getattr(args, "mu", 0.0)
    sigma = getattr(args, "sigma", 0.0)
    out_suffix = getattr(args, "out_suffix", "")

    tokens = [init_cond, geometry]
    if pdil:
        tokens.append(f"pdil{pdil:.3g}")
    if edge_weight != "flip":
        tokens.append(f"mu{mu:.3g}_s{sigma:.3g}")
    if out_suffix:
        tokens.append(out_suffix)
    return join_non_empty("_", *tokens)


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


def prepare_lattice(args: Any, **kwargs: Any) -> Any:
    """Return a cached or freshly generated Lattice3D instance (NX or GT).

    Supports ``--graph_engine`` argument: when ``args.graph_engine == 'gt'``,
    creates a ``Lattice3DGT`` instance instead of the default NX class.

    When ``init_cond`` is ``"rand"`` (or not a ground-state variant),
    skips the expensive eigenvector/energy precomputation so that large
    lattices (L>=20) can be created without full diagonalisation.
    """
    engine = get_graph_engine(args)

    if engine == "gt":
        Cls = resolve_graph_class("Lattice3D", "gt")
        dim = _normalize_dimension(getattr(args, "L"))
        l3d_kwargs = dict(
            dim=dim,
            geo=getattr(args, "geometry", None),
            pflip=getattr(args, "p", None),
            seed=getattr(args, "seed", None),
        )
        l3d_kwargs = {k: v for k, v in l3d_kwargs.items() if v is not None}
        lattice = Cls(**l3d_kwargs, **kwargs)
        lattice.flip_random_fract_edges()
        return lattice

    # Determine whether eigenspace computation is needed
    ic = getattr(args, "init_cond", "rand")
    needs_eigV = ic.startswith("ground_state") or ic.startswith("gs")

    # Default NX path (preserves load_or_compute and cell_type logic)
    base_kwargs = initialize_l3d_dict_args(args)
    cell_type = getattr(args, "cell_type", None)
    if cell_type is not None:
        base_kwargs["cell_type"] = cell_type

    if needs_eigV:
        return load_or_compute_Lattice3D(**base_kwargs, **kwargs)

    # Fast path: construct lattice without eigenspace computation
    base_kwargs.pop("cell_type", None)  # only used by load_or_compute
    lattice = Lattice3D(**base_kwargs, **kwargs)
    pflip = getattr(lattice, "pflip", 0.0)
    if pflip and pflip > 0.0:
        lattice.flip_random_fract_edges()
    return lattice
