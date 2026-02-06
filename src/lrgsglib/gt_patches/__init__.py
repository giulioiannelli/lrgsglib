"""
gt_patches: graph-tool backend for lrgsglib.

This module provides graph-tool implementations of signed graph types,
offering significant performance improvements for large-scale graphs
through graph-tool's C++ backend.

Architecture
------------
- SignedGraphGT: Base class mirroring nx_patches.SignedGraph
- Converters: Bidirectional nx <-> gt conversion utilities
- Graph types: GT implementations of lattice, random, etc.

Performance Benefits
--------------------
- graph-tool uses C++ internally with OpenMP parallelization
- Spectral operations via ARPACK/scipy are faster
- Property maps provide efficient attribute storage
- Better memory efficiency for large graphs (>100k nodes)

Example
-------
>>> from lrgsglib.gt_patches import SignedGraphGT
>>> from lrgsglib.gt_patches.converters import nx_to_gt, gt_to_nx
>>>
>>> # Convert existing nx graph
>>> gt_graph = nx_to_gt(nx_signed_graph)
>>>
>>> # Or use GT-native types
>>> from lrgsglib.gt_patches.lattice import Lattice2DGT
>>> lat = Lattice2DGT(side=100, geo='sqr')

Notes
-----
Requires graph-tool to be installed (conda install -c conda-forge graph-tool).
"""

# Check for graph-tool availability
try:
    import graph_tool
    GT_AVAILABLE = True
except ImportError:
    GT_AVAILABLE = False

if GT_AVAILABLE:
    from .signed_graph_gt import SignedGraphGT
    from .converters import nx_to_gt, gt_to_nx, GTNXConverter
    from .WeightedGraph import WeightedGraph
    from .LatticeND import LatticeND
    from .Lattice2DGT import Lattice2DGT
    from .Lattice3DGT import Lattice3DGT

    # C++ extension functions (may fail if not built)
    try:
        from .cpp import create_triangular_lattice
        _CPP_EXT_AVAILABLE = True
    except ImportError:
        _CPP_EXT_AVAILABLE = False

        def create_triangular_lattice(*args, **kwargs):
            raise ImportError(
                "C++ extensions not built. Run 'make cpp-make' from lrgsglib root."
            )

    __all__ = [
        "SignedGraphGT",
        "Lattice2DGT",
        "Lattice3DGT",
        "WeightedGraph",
        "LatticeND",
        "nx_to_gt",
        "gt_to_nx",
        "GTNXConverter",
        "GT_AVAILABLE",
        "create_triangular_lattice",
    ]
else:
    # Provide stubs that raise ImportError with helpful message
    def _gt_not_available(*args, **kwargs):
        raise ImportError(
            "graph-tool is not installed. Install with: "
            "conda install -c conda-forge graph-tool"
        )

    SignedGraphGT = _gt_not_available
    nx_to_gt = _gt_not_available
    gt_to_nx = _gt_not_available

    __all__ = ["GT_AVAILABLE"]
