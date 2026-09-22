"""
NetworkX implementations of signed graph types.

This module contains all NetworkX-based implementations with the NX suffix,
plus utility modules (funcs, datasets) and supporting functions.

Usage
-----
>>> from lrgsglib.graphs.nx import Lattice2DNX, ErdosRenyiNX
>>> from lrgsglib.graphs.nx.funcs import signed_laplacian_matrix
>>> from lrgsglib.graphs.nx.datasets import RealDatasetLoader

Categories
----------
Base:
    SignedGraphNX : Base signed graph class

Lattice:
    LatticeNDNX : N-dimensional cubic lattice
    Lattice2DNX : 2D lattice graph
    Lattice3DNX : 3D lattice graph

Random:
    RandomGraphNX, ErdosRenyiNX, BarabasiAlbertNX, WattsStrogatzNX, etc.

Complete:
    CompleteGraphNX, FullyConnectedNX

Neural:
    HofieldNNNX, SCSGeneralizedNNNX

Fractal:
    FractalGraphNX, DGMgraphNX, SierpinskiNX

Bipartite:
    BipartiteGraphNX (unified: random + degree-sequence modes), BipartiteFromDegreeSequenceNX (alias)

Multispectral:
    MultispectralGraphNX, MultiplicativeCascadeGraphNX, VicsekGraphNX, HierarchicalModularNetworkNX

Dirac:
    DiracLatticeGraphNX, DiracCombGraphNX, DiracBrushGraphNX

Temporal:
    TemporalGraphNX, TemporalSignedGraphNX
"""

# Base
from .SignedGraphNX import SignedGraphNX

SignedGraph = SignedGraphNX

from .Lattice2DNX import (
    Lattice2D,
    Lattice2DNX,
    create_lattice_with_eigenspace,
    load_or_compute_Lattice2D,
    load_or_compute_Lattice2DNX,
)
from .Lattice3DNX import (
    Lattice3D,
    Lattice3DNX,
    load_or_compute_Lattice3D,
    load_or_compute_Lattice3DNX,
)

# Lattice
from .LatticeNDNX import LatticeNDNX

# Random
from .RandomGraphNX import RandomGraphNX

RandomGraph = RandomGraphNX
from .BarabasiAlbertNX import BarabasiAlbert, BarabasiAlbertNX
from .ErdosRenyiNX import ErdosRenyi, ErdosRenyiNX
from .WattsStrogatzNX import WattsStrogatz, WattsStrogatzNX

WattStrogatz = WattsStrogatzNX  # legacy typo alias
from .BipartiteFromDegreeSequenceNX import (
    BipartiteFromDegreeSequence,
    BipartiteFromDegreeSequenceNX,
)

# Bipartite
from .BipartiteGraphNX import BipartiteGraph, BipartiteGraphNX

# Complete
from .CompleteGraphNX import CompleteGraph, CompleteGraphNX
from .ConfigurationModelNX import ConfigurationModel, ConfigurationModelNX
from .DGMgraphNX import DGMgraph, DGMgraphNX

# Dirac
from .DiracLatticeNX import (
    DiracBrushGraph,
    DiracBrushGraphNX,
    DiracCombGraph,
    DiracCombGraphNX,
    DiracLatticeGraph,
    DiracLatticeGraphNX,
)
from .DualBarabasiAlbertNX import DualBarabasiAlbert, DualBarabasiAlbertNX
from .ExtendedBarabasiAlbertNX import (
    ExtendedBarabasiAlbert,
    ExtendedBarabasiAlbertNX,
)

# Fractal
from .FractalGraphNX import FractalGraph, FractalGraphNX
from .FullyConnectedNX import FullyConnected, FullyConnectedNX
from .HierarchicalModularNX import (
    HierarchicalModularNetwork,
    HierarchicalModularNetworkNX,
)

# Neural
from .HofieldNNNX import HofieldNN, HofieldNNNX
from .HolmeKimNX import HolmeKim, HolmeKimNX
from .kRegularGraphNX import kRegularGraph, kRegularGraphNX
from .LFRBenchmarkNX import LFRBenchmark, LFRBenchmarkNX
from .MultiplicativeCascadeNX import (
    MultiplicativeCascadeGraph,
    MultiplicativeCascadeGraphNX,
)

# Multispectral
from .MultispectralGraphNX import MultispectralGraph, MultispectralGraphNX
from .RandomGeometricNX import RandomGeometric, RandomGeometricNX
from .SCSGeneralizedNNNX import SCSGeneralizedNN, SCSGeneralizedNNNX
from .SierpinskiNX import SierpinskiGraph, SierpinskiGraphNX, SierpinskiNX
from .StochasticBlockModelNX import StochasticBlockModel, StochasticBlockModelNX
from .VicsekNX import VicsekGraph, VicsekGraphNX

DiracLattice = DiracLatticeGraphNX

# Utility submodules
from . import datasets, funcs
from .datasets import BUILTIN_DATASETS, RealDatasetLoader

# GraphOfGraphs
from .GraphOfGraphsNX import GraphOfGraphs as GraphOfGraphsNXAlias
from .GraphOfGraphsNX import GraphOfGraphsNX

# Temporal
from .TemporalGraphNX import TemporalGraph, TemporalGraphNX
from .TemporalSignedGraphNX import TemporalSignedGraph, TemporalSignedGraphNX

__all__ = [
    # Base
    "SignedGraphNX",
    "SignedGraph",
    # Lattice
    "LatticeNDNX",
    "Lattice2DNX",
    "Lattice2D",
    "Lattice3DNX",
    "Lattice3D",
    "create_lattice_with_eigenspace",
    "load_or_compute_Lattice2DNX",
    "load_or_compute_Lattice2D",
    "load_or_compute_Lattice3DNX",
    "load_or_compute_Lattice3D",
    # Random
    "RandomGraphNX",
    "RandomGraph",
    "ErdosRenyiNX",
    "ErdosRenyi",
    "BarabasiAlbertNX",
    "BarabasiAlbert",
    "WattsStrogatzNX",
    "WattsStrogatz",
    "WattStrogatz",
    "StochasticBlockModelNX",
    "StochasticBlockModel",
    "kRegularGraphNX",
    "kRegularGraph",
    "ConfigurationModelNX",
    "ConfigurationModel",
    "RandomGeometricNX",
    "RandomGeometric",
    "LFRBenchmarkNX",
    "LFRBenchmark",
    "ExtendedBarabasiAlbertNX",
    "ExtendedBarabasiAlbert",
    "DualBarabasiAlbertNX",
    "DualBarabasiAlbert",
    "HolmeKimNX",
    "HolmeKim",
    # Complete
    "CompleteGraphNX",
    "CompleteGraph",
    "FullyConnectedNX",
    "FullyConnected",
    # Neural
    "HofieldNNNX",
    "HofieldNN",
    "SCSGeneralizedNNNX",
    "SCSGeneralizedNN",
    # Fractal
    "FractalGraphNX",
    "FractalGraph",
    "DGMgraphNX",
    "DGMgraph",
    "SierpinskiNX",
    "SierpinskiGraph",
    # Bipartite
    "BipartiteGraphNX",
    "BipartiteGraph",
    "BipartiteFromDegreeSequenceNX",
    "BipartiteFromDegreeSequence",
    # Multispectral
    "MultispectralGraphNX",
    "MultispectralGraph",
    "MultiplicativeCascadeGraphNX",
    "MultiplicativeCascadeGraph",
    "VicsekGraphNX",
    "VicsekGraph",
    "HierarchicalModularNetworkNX",
    "HierarchicalModularNetwork",
    # Dirac
    "DiracLatticeGraphNX",
    "DiracLatticeGraph",
    "DiracBrushGraphNX",
    "DiracBrushGraph",
    "DiracCombGraphNX",
    "DiracCombGraph",
    "DiracLattice",
    # GraphOfGraphs
    "GraphOfGraphsNX",
    # Temporal
    "TemporalGraphNX",
    "TemporalGraph",
    "TemporalSignedGraphNX",
    "TemporalSignedGraph",
    # Datasets
    "RealDatasetLoader",
    "BUILTIN_DATASETS",
    # Submodules
    "funcs",
    "datasets",
]
