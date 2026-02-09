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
    LatticeNDNX : N-dimensional lattice (abstract base)
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

# Lattice
from .Lattice2DNX import (
    Lattice2DNX,
    Lattice2D,
    create_lattice_with_eigenspace,
    load_or_compute_Lattice2DNX,
    load_or_compute_Lattice2D,
)
from .Lattice3DNX import (
    Lattice3DNX,
    Lattice3D,
    load_or_compute_Lattice3DNX,
    load_or_compute_Lattice3D,
)

# Random
from .RandomGraphNX import RandomGraphNX
RandomGraph = RandomGraphNX
from .ErdosRenyiNX import ErdosRenyiNX, ErdosRenyi
from .BarabasiAlbertNX import BarabasiAlbertNX, BarabasiAlbert
from .WattsStrogatzNX import WattsStrogatzNX, WattsStrogatz
WattStrogatz = WattsStrogatzNX  # legacy typo alias
from .StochasticBlockModelNX import StochasticBlockModelNX, StochasticBlockModel
from .kRegularGraphNX import kRegularGraphNX, kRegularGraph
from .ConfigurationModelNX import ConfigurationModelNX, ConfigurationModel
from .RandomGeometricNX import RandomGeometricNX, RandomGeometric
from .LFRBenchmarkNX import LFRBenchmarkNX, LFRBenchmark
from .ExtendedBarabasiAlbertNX import ExtendedBarabasiAlbertNX, ExtendedBarabasiAlbert
from .DualBarabasiAlbertNX import DualBarabasiAlbertNX, DualBarabasiAlbert
from .HolmeKimNX import HolmeKimNX, HolmeKim

# Complete
from .CompleteGraphNX import CompleteGraphNX, CompleteGraph
from .FullyConnectedNX import FullyConnectedNX, FullyConnected

# Neural
from .HofieldNNNX import HofieldNNNX, HofieldNN
from .SCSGeneralizedNNNX import SCSGeneralizedNNNX, SCSGeneralizedNN

# Fractal
from .FractalGraphNX import FractalGraphNX, FractalGraph
from .DGMgraphNX import DGMgraphNX, DGMgraph
from .SierpinskiNX import SierpinskiNX, SierpinskiGraph, SierpinskiGraphNX

# Bipartite
from .BipartiteGraphNX import BipartiteGraphNX, BipartiteGraph
from .BipartiteFromDegreeSequenceNX import BipartiteFromDegreeSequenceNX, BipartiteFromDegreeSequence

# Multispectral
from .MultispectralGraphNX import MultispectralGraphNX, MultispectralGraph
from .MultiplicativeCascadeNX import MultiplicativeCascadeGraphNX, MultiplicativeCascadeGraph
from .VicsekNX import VicsekGraphNX, VicsekGraph
from .HierarchicalModularNX import HierarchicalModularNetworkNX, HierarchicalModularNetwork

# Dirac
from .DiracLatticeNX import (
    DiracLatticeGraphNX,
    DiracLatticeGraph,
    DiracBrushGraphNX,
    DiracBrushGraph,
    DiracCombGraphNX,
    DiracCombGraph,
)
DiracLattice = DiracLatticeGraphNX

# GraphOfGraphs
from .GraphOfGraphsNX import (
    GraphOfGraphsNX,
    GraphOfGraphs as GraphOfGraphsNXAlias,
)

# Temporal
from .TemporalGraphNX import TemporalGraphNX, TemporalGraph
from .TemporalSignedGraphNX import TemporalSignedGraphNX, TemporalSignedGraph

# Utility submodules
from . import funcs
from . import datasets
from .datasets import RealDatasetLoader, BUILTIN_DATASETS

__all__ = [
    # Base
    "SignedGraphNX",
    "SignedGraph",
    # Lattice
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
