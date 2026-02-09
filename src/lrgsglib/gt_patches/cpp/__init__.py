"""
gt_patches.cpp — DEPRECATED.

C++ extensions have moved to per-class ``cpp/`` subdirectories under
``lrgsglib.graphs.gt``.  E.g. triangular lattice is now at
``lrgsglib.graphs.gt.Lattice2DGT.cpp``.
"""
import warnings as _warnings

_WARNED = False


def __getattr__(name):
    global _WARNED
    if name == "create_triangular_lattice":
        if not _WARNED:
            _warnings.warn(
                "Importing 'create_triangular_lattice' from "
                "'lrgsglib.gt_patches.cpp' is deprecated. "
                "Use 'lrgsglib.graphs.gt.Lattice2DGT.cpp' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            _WARNED = True
        from lrgsglib.graphs.gt.Lattice2DGT.cpp import (
            create_triangular_lattice,
        )
        return create_triangular_lattice

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__: list[str] = []
