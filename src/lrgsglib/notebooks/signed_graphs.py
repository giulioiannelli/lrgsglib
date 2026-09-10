"""Signed-graph constructors and cached lattice loaders for notebooks."""
# Engine-agnostic graph factories (preferred interface in notebooks).
from ..graphs import (
    SignedGraph,
    Disorder,
    CompositeDisorder,
    register_coupling,
    register_support,
    registered_supports,
    Lattice2D,
    Lattice3D,
    LatticeND,
    ErdosRenyi,
    StochasticBlockModel,
    HierarchicalModular,
    MultiplicativeCascade,
    VicsekGraph,
    DiracCombGraph,
    DiracBrushGraph,
    SierpinskiGraph,
    BarabasiAlbert,
    WattsStrogatz,
    FullyConnected,
    kRegularGraph,
    BipartiteGraph,
    RandomGeometric,
    ConfigurationModel,
    LFRBenchmark,
    HolmeKim,
    DualBarabasiAlbert,
    ExtendedBarabasiAlbert,
    DGMgraph,
)

# Cached-eigenspace lattice loaders (NX engine), part of the historical
# notebook API surface (same names as the top-level lrgsglib exports).
from ..graphs.nx import (
    load_or_compute_Lattice2D,
    load_or_compute_Lattice3D,
)
