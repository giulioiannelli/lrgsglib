"""
graph-tool implementations of signed graph types.

This module contains all graph-tool-based implementations with the GT suffix.

Categories
----------
Base:
    SignedGraphGT : Base signed graph class (graph-tool backend)

Lattice:
    LatticeNDGT : N-dimensional lattice (native GT)
    Lattice2DGT : 2D lattice graph (native GT)
    Lattice3DGT : 3D lattice graph (native GT)

Random:
    RandomGraphGT : Random graph base class (fallback to NX)
    ErdosRenyiGT : Erdos-Renyi random graph (native GT)
    BarabasiAlbertGT : Barabasi-Albert scale-free graph (native GT)
    WattsStrogatzGT : Watts-Strogatz small-world graph (native GT)
    StochasticBlockModelGT : Stochastic Block Model graph (native GT)
    kRegularGraphGT : k-regular random graph (fallback to NX)
    ConfigurationModelGT : Configuration model graph (fallback to NX)
    RandomGeometricGT : Random geometric graph (fallback to NX)
    LFRBenchmarkGT : LFR benchmark graph (fallback to NX)
    ExtendedBarabasiAlbertGT : Extended BA (fallback to NX)
    DualBarabasiAlbertGT : Dual BA (fallback to NX)
    HolmeKimGT : BA with clustering (fallback to NX)

Complete:
    CompleteGraphGT : Complete graph (fallback to NX)
    FullyConnectedGT : Fully connected graph (fallback to NX)

Neural:
    HofieldNNGT : Hopfield neural network graph (fallback to NX)
    SCSGeneralizedNNGT : SCS generalized neural network graph (fallback to NX)

Fractal:
    FractalGraphGT : Fractal graph base class (fallback to NX)
    DGMgraphGT : Dorogovtsev-Goltsev-Mendes graph (fallback to NX)
    SierpinskiGraphGT : Sierpinski fractal graph (fallback to NX)

Bipartite:
    BipartiteGraphGT : Bipartite graph (fallback to NX)
    BipartiteFromDegreeSequenceGT : Bipartite from degree sequence (fallback to NX)

Multispectral:
    MultispectralGraphGT : Multispectral graph (fallback to NX)
    MultiplicativeCascadeGraphGT : Multiplicative cascade graph (fallback to NX)
    VicsekGraphGT : Vicsek fractal graph (fallback to NX)
    HierarchicalModularNetworkGT : Hierarchical modular network (fallback to NX)

Dirac:
    DiracLatticeGraphGT : Dirac lattice graph (fallback to NX)
    DiracCombGraphGT : Dirac comb graph (fallback to NX)
    DiracBrushGraphGT : Dirac brush graph (fallback to NX)

Temporal:
    TemporalGraphGT : Temporal graph (fallback to NX)
    TemporalSignedGraphGT : Temporal signed graph (fallback to NX)

Note: Classes marked "(fallback to NX)" will emit a warning when used and
delegate to the NetworkX implementation until native GT implementations are added.
"""

# Lazy imports to handle missing graph-tool gracefully
def __getattr__(name):
    # Base
    if name == "SignedGraphGT":
        from .base import SignedGraphGT
        return SignedGraphGT

    # Lattice (native GT)
    elif name == "LatticeNDGT":
        from .lattice import LatticeNDGT
        return LatticeNDGT
    elif name == "Lattice2DGT":
        from .lattice import Lattice2DGT
        return Lattice2DGT
    elif name == "Lattice3DGT":
        from .lattice import Lattice3DGT
        return Lattice3DGT

    # Random (native GT)
    elif name == "ErdosRenyiGT":
        from .random import ErdosRenyiGT
        return ErdosRenyiGT
    elif name == "BarabasiAlbertGT":
        from .random import BarabasiAlbertGT
        return BarabasiAlbertGT
    elif name == "WattsStrogatzGT":
        from .random import WattsStrogatzGT
        return WattsStrogatzGT
    elif name == "StochasticBlockModelGT":
        from .random import StochasticBlockModelGT
        return StochasticBlockModelGT

    # Random (fallback to NX)
    elif name == "RandomGraphGT":
        from .random import RandomGraphGT
        return RandomGraphGT
    elif name == "kRegularGraphGT":
        from .random import kRegularGraphGT
        return kRegularGraphGT
    elif name == "ConfigurationModelGT":
        from .random import ConfigurationModelGT
        return ConfigurationModelGT
    elif name == "RandomGeometricGT":
        from .random import RandomGeometricGT
        return RandomGeometricGT
    elif name == "LFRBenchmarkGT":
        from .random import LFRBenchmarkGT
        return LFRBenchmarkGT
    elif name == "ExtendedBarabasiAlbertGT":
        from .random import ExtendedBarabasiAlbertGT
        return ExtendedBarabasiAlbertGT
    elif name == "DualBarabasiAlbertGT":
        from .random import DualBarabasiAlbertGT
        return DualBarabasiAlbertGT
    elif name == "HolmeKimGT":
        from .random import HolmeKimGT
        return HolmeKimGT

    # Complete (fallback to NX)
    elif name == "CompleteGraphGT":
        from .complete import CompleteGraphGT
        return CompleteGraphGT
    elif name == "FullyConnectedGT":
        from .complete import FullyConnectedGT
        return FullyConnectedGT

    # Neural (fallback to NX)
    elif name == "HofieldNNGT":
        from .neural import HofieldNNGT
        return HofieldNNGT
    elif name == "SCSGeneralizedNNGT":
        from .neural import SCSGeneralizedNNGT
        return SCSGeneralizedNNGT

    # Fractal (fallback to NX)
    elif name == "FractalGraphGT":
        from .fractal import FractalGraphGT
        return FractalGraphGT
    elif name == "DGMgraphGT":
        from .fractal import DGMgraphGT
        return DGMgraphGT
    elif name == "SierpinskiGraphGT":
        from .fractal import SierpinskiGraphGT
        return SierpinskiGraphGT

    # Bipartite (fallback to NX)
    elif name == "BipartiteGraphGT":
        from .bipartite import BipartiteGraphGT
        return BipartiteGraphGT
    elif name == "BipartiteFromDegreeSequenceGT":
        from .bipartite import BipartiteFromDegreeSequenceGT
        return BipartiteFromDegreeSequenceGT

    # Multispectral (fallback to NX)
    elif name == "MultispectralGraphGT":
        from .multispectral import MultispectralGraphGT
        return MultispectralGraphGT
    elif name == "MultiplicativeCascadeGraphGT":
        from .multispectral import MultiplicativeCascadeGraphGT
        return MultiplicativeCascadeGraphGT
    elif name == "VicsekGraphGT":
        from .multispectral import VicsekGraphGT
        return VicsekGraphGT
    elif name == "HierarchicalModularNetworkGT":
        from .multispectral import HierarchicalModularNetworkGT
        return HierarchicalModularNetworkGT

    # Dirac (fallback to NX)
    elif name == "DiracLatticeGraphGT":
        from .dirac import DiracLatticeGraphGT
        return DiracLatticeGraphGT
    elif name == "DiracCombGraphGT":
        from .dirac import DiracCombGraphGT
        return DiracCombGraphGT
    elif name == "DiracBrushGraphGT":
        from .dirac import DiracBrushGraphGT
        return DiracBrushGraphGT

    # Temporal (fallback to NX)
    elif name == "TemporalGraphGT":
        from .temporal import TemporalGraphGT
        return TemporalGraphGT
    elif name == "TemporalSignedGraphGT":
        from .temporal import TemporalSignedGraphGT
        return TemporalSignedGraphGT

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    # Base
    "SignedGraphGT",
    # Lattice (native GT)
    "LatticeNDGT",
    "Lattice2DGT",
    "Lattice3DGT",
    # Random (native GT)
    "ErdosRenyiGT",
    "BarabasiAlbertGT",
    "WattsStrogatzGT",
    "StochasticBlockModelGT",
    # Random (fallback to NX)
    "RandomGraphGT",
    "kRegularGraphGT",
    "ConfigurationModelGT",
    "RandomGeometricGT",
    "LFRBenchmarkGT",
    "ExtendedBarabasiAlbertGT",
    "DualBarabasiAlbertGT",
    "HolmeKimGT",
    # Complete (fallback to NX)
    "CompleteGraphGT",
    "FullyConnectedGT",
    # Neural (fallback to NX)
    "HofieldNNGT",
    "SCSGeneralizedNNGT",
    # Fractal (fallback to NX)
    "FractalGraphGT",
    "DGMgraphGT",
    "SierpinskiGraphGT",
    # Bipartite (fallback to NX)
    "BipartiteGraphGT",
    "BipartiteFromDegreeSequenceGT",
    # Multispectral (fallback to NX)
    "MultispectralGraphGT",
    "MultiplicativeCascadeGraphGT",
    "VicsekGraphGT",
    "HierarchicalModularNetworkGT",
    # Dirac (fallback to NX)
    "DiracLatticeGraphGT",
    "DiracCombGraphGT",
    "DiracBrushGraphGT",
    # Temporal (fallback to NX)
    "TemporalGraphGT",
    "TemporalSignedGraphGT",
]
