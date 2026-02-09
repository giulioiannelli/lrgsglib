"""
gt_patches — DEPRECATED.

All graph-tool implementations have moved to ``lrgsglib.graphs.gt``.
This shim re-exports public names for backward compatibility but will
emit a DeprecationWarning on first access.

Migration guide::

    # Old
    from lrgsglib.gt_patches import SignedGraphGT
    from lrgsglib.gt_patches.converters import nx_to_gt

    # New
    from lrgsglib.graphs.gt import SignedGraphGT
    from lrgsglib.graphs.gt import nx_to_gt
"""
import warnings as _warnings

# Check for graph-tool availability (kept for backward compat)
try:
    import graph_tool
    GT_AVAILABLE = True
except ImportError:
    GT_AVAILABLE = False

_WARNED = set()


def __getattr__(name):
    # Only shim names that used to live here
    _PUBLIC = {
        "SignedGraphGT", "Lattice2DGT", "Lattice3DGT",
        "LatticeND", "WeightedGraph",
        "nx_to_gt", "gt_to_nx", "GTNXConverter",
        "create_triangular_lattice",
    }
    if name in _PUBLIC:
        if name not in _WARNED:
            _warnings.warn(
                f"Importing {name!r} from 'lrgsglib.gt_patches' is deprecated. "
                f"Use 'lrgsglib.graphs.gt' instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            _WARNED.add(name)

        # Resolve from new canonical location
        from lrgsglib.graphs.gt import _converters

        if name == "SignedGraphGT":
            from lrgsglib.graphs.gt.SignedGraphGT import SignedGraphGT
            return SignedGraphGT
        elif name == "Lattice2DGT":
            from lrgsglib.graphs.gt.Lattice2DGT import Lattice2DGT
            return Lattice2DGT
        elif name == "Lattice3DGT":
            from lrgsglib.graphs.gt.Lattice3DGT import Lattice3DGT
            return Lattice3DGT
        elif name == "LatticeND":
            from lrgsglib.graphs.gt.LatticeNDGT import LatticeNDGT
            return LatticeNDGT
        elif name == "WeightedGraph":
            # WeightedGraph removed — return Graph directly
            if not GT_AVAILABLE:
                raise ImportError("graph-tool is not installed")
            return graph_tool.Graph
        elif name in ("nx_to_gt", "gt_to_nx", "GTNXConverter"):
            return getattr(_converters, name)
        elif name == "create_triangular_lattice":
            from lrgsglib.graphs.gt.Lattice2DGT.cpp import (
                create_triangular_lattice,
            )
            return create_triangular_lattice

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ["GT_AVAILABLE"]
