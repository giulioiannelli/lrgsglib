"""Signed-graph constructors and cached lattice loaders for notebooks."""

# Engine-agnostic graph factories (preferred interface in notebooks).
from ..graphs import (
    BarabasiAlbert,
    BipartiteGraph,
    CompositeDisorder,
    ConfigurationModel,
    DGMgraph,
    DiracBrushGraph,
    DiracCombGraph,
    Disorder,
    DualBarabasiAlbert,
    ErdosRenyi,
    ExtendedBarabasiAlbert,
    FullyConnected,
    HierarchicalModular,
    HolmeKim,
    Lattice2D,
    Lattice3D,
    LatticeND,
    LFRBenchmark,
    MultiplicativeCascade,
    RandomGeometric,
    SierpinskiGraph,
    SignedGraph,
    StochasticBlockModel,
    VicsekGraph,
    WattsStrogatz,
    kRegularGraph,
    register_coupling,
    register_support,
    registered_supports,
)

# Cached-eigenspace lattice loaders (NX engine), part of the historical
# notebook API surface (same names as the top-level lrgsglib exports).
from ..graphs.nx import (
    load_or_compute_Lattice2D,
    load_or_compute_Lattice3D,
)
