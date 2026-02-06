"""
DEPRECATED: nx_patches module.

This module is deprecated. Please use lrgsglib.graphs.nx instead.

All imports from this module will continue to work but will emit deprecation
warnings. The implementations have been moved to graphs/nx/ with the NX suffix
naming convention.

Migration Examples
------------------
Old (deprecated):
    >>> from lrgsglib.nx_patches import Lattice2D, ErdosRenyi, SignedGraph

New (recommended):
    >>> from lrgsglib.graphs.nx import Lattice2DNX, ErdosRenyiNX, SignedGraphNX
    # Or use backward compatibility aliases:
    >>> from lrgsglib.graphs.nx import Lattice2D, ErdosRenyi

The unified graphs module provides both NX and GT backends:
    >>> from lrgsglib.graphs import Lattice2D
    >>> lat = Lattice2D(10, engine='nx')  # NetworkX backend
    >>> lat = Lattice2D(10, engine='gt')  # graph-tool backend
"""

import warnings
import sys

# Emit deprecation warning on import
warnings.warn(
    "Importing from lrgsglib.nx_patches is deprecated. "
    "Use lrgsglib.graphs.nx instead. "
    "See lrgsglib.graphs for the unified multi-engine API.",
    DeprecationWarning,
    stacklevel=2
)

# Re-export everything from graphs.nx with original names for backward compatibility
from ..graphs.nx import (
    # Base
    SignedGraphNX as SignedGraph,
    # Lattice
    Lattice2DNX as Lattice2D,
    Lattice3DNX as Lattice3D,
    load_or_compute_Lattice2DNX as load_or_compute_Lattice2D,
    load_or_compute_Lattice3DNX as load_or_compute_Lattice3D,
    create_lattice_with_eigenspace,
    # Random
    RandomGraphNX as RandomGraph,
    ErdosRenyiNX as ErdosRenyi,
    BarabasiAlbertNX as BarabasiAlbert,
    WattsStrogatzNX as WattsStrogatz,
    WattsStrogatzNX as WattStrogatz,  # typo alias for legacy code
    StochasticBlockModelNX as StochasticBlockModel,
    kRegularGraphNX as kRegularGraph,
    ConfigurationModelNX as ConfigurationModel,
    RandomGeometricNX as RandomGeometric,
    LFRBenchmarkNX as LFRBenchmark,
    ExtendedBarabasiAlbertNX as ExtendedBarabasiAlbert,
    DualBarabasiAlbertNX as DualBarabasiAlbert,
    HolmeKimNX as HolmeKim,
    # Complete
    CompleteGraphNX as CompleteGraph,
    FullyConnectedNX as FullyConnected,
    # Neural
    HofieldNNNX as HofieldNN,
    SCSGeneralizedNNNX as SCSGeneralizedNN,
    # Fractal
    FractalGraphNX as FractalGraph,
    DGMgraphNX as DGMgraph,
    SierpinskiNX as SierpinskiGraph,
    # Bipartite
    BipartiteGraphNX as BipartiteGraph,
    BipartiteFromDegreeSequenceNX as BipartiteFromDegreeSequence,
    # Multispectral
    MultispectralGraphNX as MultispectralGraph,
    MultiplicativeCascadeGraphNX as MultiplicativeCascadeGraph,
    VicsekGraphNX as VicsekGraph,
    HierarchicalModularNetworkNX as HierarchicalModularNetwork,
    # Dirac
    DiracLatticeGraphNX as DiracLatticeGraph,
    DiracLatticeGraphNX as DiracLattice,  # alias
    DiracBrushGraphNX as DiracBrushGraph,
    DiracCombGraphNX as DiracCombGraph,
    # Temporal
    TemporalGraphNX as TemporalGraph,
    TemporalSignedGraphNX as TemporalSignedGraph,
    # Utility submodules
    funcs,
    datasets,
)

# Re-export dataset utilities
from ..graphs.nx.datasets import RealDatasetLoader, BUILTIN_DATASETS

# LatticeND needs special handling - check if it exists
try:
    from ..graphs.nx import LatticeNDNX as LatticeND
except ImportError:
    LatticeND = None

__all__ = [
    # Base
    "SignedGraph",
    # Lattice
    "Lattice2D",
    "Lattice3D",
    "LatticeND",
    "load_or_compute_Lattice2D",
    "load_or_compute_Lattice3D",
    "create_lattice_with_eigenspace",
    # Random
    "RandomGraph",
    "ErdosRenyi",
    "BarabasiAlbert",
    "WattsStrogatz",
    "WattStrogatz",
    "StochasticBlockModel",
    "kRegularGraph",
    "ConfigurationModel",
    "RandomGeometric",
    "LFRBenchmark",
    "ExtendedBarabasiAlbert",
    "DualBarabasiAlbert",
    "HolmeKim",
    # Complete
    "CompleteGraph",
    "FullyConnected",
    # Neural
    "HofieldNN",
    "SCSGeneralizedNN",
    # Fractal
    "FractalGraph",
    "DGMgraph",
    "SierpinskiGraph",
    # Bipartite
    "BipartiteGraph",
    "BipartiteFromDegreeSequence",
    # Multispectral
    "MultispectralGraph",
    "MultiplicativeCascadeGraph",
    "VicsekGraph",
    "HierarchicalModularNetwork",
    # Dirac
    "DiracLatticeGraph",
    "DiracLattice",
    "DiracBrushGraph",
    "DiracCombGraph",
    # Temporal
    "TemporalGraph",
    "TemporalSignedGraph",
    # Datasets
    "RealDatasetLoader",
    "BUILTIN_DATASETS",
    # Submodules
    "funcs",
    "datasets",
]

# Store class references to override submodule shadowing
# When Python imports a submodule like `lrgsglib.nx_patches.SignedGraph._backend`,
# it caches `SignedGraph` as a submodule (folder), which shadows the class in __init__.py.
# This mapping ensures the correct class is returned via custom __getattribute__.
_CLASS_OVERRIDES = {
    "SignedGraph": SignedGraph,
    "Lattice2D": Lattice2D,
    "Lattice3D": Lattice3D,
    "ErdosRenyi": ErdosRenyi,
    "FullyConnected": FullyConnected,
    "HofieldNN": HofieldNN,
    "SCSGeneralizedNN": SCSGeneralizedNN,
    "DGMgraph": DGMgraph,
    "MultispectralGraph": MultispectralGraph,
    "MultiplicativeCascadeGraph": MultiplicativeCascadeGraph,
    "VicsekGraph": VicsekGraph,
    "HierarchicalModularNetwork": HierarchicalModularNetwork,
    "DiracLatticeGraph": DiracLatticeGraph,
    "DiracLattice": DiracLattice,
    "DiracBrushGraph": DiracBrushGraph,
    "DiracCombGraph": DiracCombGraph,
}


# Use a custom module class to intercept attribute access
# This is needed because __getattr__ is only called when the attribute
# is NOT found in __dict__, but Python adds submodules to __dict__
import types


class _NXPatchesModule(types.ModuleType):
    """Custom module type that overrides submodule shadowing."""

    _CLASS_OVERRIDES_ATTR = {}  # Will be set after instantiation

    def __getattribute__(self, name: str):
        # Use object.__getattribute__ to avoid recursion
        try:
            overrides = object.__getattribute__(self, "_CLASS_OVERRIDES_ATTR")
        except AttributeError:
            overrides = {}
        if name in overrides:
            return overrides[name]
        return super().__getattribute__(name)


# Replace this module with our custom module type
_current_module = sys.modules[__name__]
_new_module = _NXPatchesModule(__name__, __doc__)
# Set overrides BEFORE updating dict to avoid issues
_new_module._CLASS_OVERRIDES_ATTR = _CLASS_OVERRIDES
_new_module.__dict__.update(_current_module.__dict__)
sys.modules[__name__] = _new_module
