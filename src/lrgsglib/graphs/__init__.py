"""
Unified graph interface with multi-engine support.

This module provides a unified API for creating signed graphs using different
backends (NetworkX, graph-tool, igraph). Users can specify the engine explicitly
or rely on a configurable global default.

Quick Start
-----------
>>> from lrgsglib.graphs import Lattice2D, set_default_engine
>>>
>>> # Use specific engine
>>> lat_nx = Lattice2D(side1=100, geo='sqr', engine='nx')
>>> lat_gt = Lattice2D(side1=100, geo='sqr', engine='gt')
>>>
>>> # Or set global default
>>> set_default_engine('gt')
>>> lat = Lattice2D(side1=100, geo='sqr')  # Uses graph-tool

Available Graph Types
--------------------
Base:
    SignedGraph : Base signed graph wrapper

Lattice:
    Lattice2D, Lattice3D, LatticeND

Random:
    ErdosRenyi, kRegularGraph, ConfigurationModel,
    RandomGeometric, LFRBenchmark

Scale-Free:
    BarabasiAlbert, ExtendedBarabasiAlbert, DualBarabasiAlbert, HolmeKim

Small-World:
    WattsStrogatz

Community:
    StochasticBlockModel

Complete:
    FullyConnected

Bipartite:
    BipartiteGraph, BipartiteFromDegreeSequence

Neural:
    HofieldNN, SCSGeneralizedNN

Multispectral:
    MultiplicativeCascade, VicsekGraph, HierarchicalModular

Fractal:
    SierpinskiGraph, DGMgraph

Graph-of-Graphs:
    GraphOfGraphs, DiracCombGraph, DiracBrushGraph

Temporal:
    TemporalGraph

Engine Configuration
--------------------
Engines can be configured via:
1. Explicit parameter: ``Lattice2D(..., engine='gt')``
2. Global setting: ``set_default_engine('gt')``
3. Environment variable: ``LRGSG_GRAPH_ENGINE=gt``

Priority: explicit parameter > environment variable > global setting

Supported Engines
----------------
- 'nx' : NetworkX (default, always available)
- 'gt' : graph-tool (high performance, requires installation)
- 'ig' : igraph (future support)
"""

from ._engine import (
    GraphEngine,
    available_engines,
    get_default_engine,
    get_implementation,
    is_engine_available,
    list_engines_for_type,
    list_registered_types,
    register_implementation,
    set_default_engine,
)

# Disorder model (engine-neutral): support x coupling-law spec carried on sg
from ._shared._disorder import (
    CompositeDisorder,
    Disorder,
    register_coupling,
    register_support,
    registered_supports,
)

# Scale-Free
from .BarabasiAlbert import BarabasiAlbert
from .BipartiteFromDegreeSequence import BipartiteFromDegreeSequence

# Bipartite
from .BipartiteGraph import BipartiteGraph
from .ConfigurationModel import ConfigurationModel
from .DGMgraph import DGMgraph
from .DiracBrushGraph import DiracBrushGraph
from .DiracCombGraph import DiracCombGraph
from .DualBarabasiAlbert import DualBarabasiAlbert

# Random
from .ErdosRenyi import ErdosRenyi
from .ExtendedBarabasiAlbert import ExtendedBarabasiAlbert

# Complete
from .FullyConnected import FullyConnected

# Graph-of-Graphs / Dirac
from .GraphOfGraphs import GraphOfGraphs
from .HierarchicalModular import HierarchicalModular

# Neural
from .HofieldNN import HofieldNN
from .HolmeKim import HolmeKim
from .kRegularGraph import kRegularGraph

# Lattice
from .Lattice2D import Lattice2D
from .Lattice3D import Lattice3D
from .LatticeND import LatticeND
from .LFRBenchmark import LFRBenchmark

# Multispectral (concrete classes)
from .MultiplicativeCascade import MultiplicativeCascade
from .protocols import (
    LatticeGraphProtocol,
    SignedGraphProtocol,
    SpectralGraphProtocol,
    is_lattice_graph,
    is_signed_graph,
)
from .RandomGeometric import RandomGeometric
from .SCSGeneralizedNN import SCSGeneralizedNN

# Fractal (concrete classes)
from .SierpinskiGraph import SierpinskiGraph

# Base
from .SignedGraph import SignedGraph

# Community
from .StochasticBlockModel import StochasticBlockModel

# Temporal
from .TemporalGraph import TemporalGraph
from .VicsekGraph import VicsekGraph

# Small-World
from .WattsStrogatz import WattsStrogatz

# === Graph type facades ===


__all__ = [
    # Engine management
    "GraphEngine",
    "set_default_engine",
    "get_default_engine",
    "get_implementation",
    "register_implementation",
    "is_engine_available",
    "available_engines",
    "list_registered_types",
    "list_engines_for_type",
    # Protocols
    "SignedGraphProtocol",
    "SpectralGraphProtocol",
    "LatticeGraphProtocol",
    "is_signed_graph",
    "is_lattice_graph",
    # Base
    "SignedGraph",
    # Disorder model
    "Disorder",
    "CompositeDisorder",
    "register_coupling",
    "register_support",
    "registered_supports",
    # Lattice
    "Lattice2D",
    "Lattice3D",
    "LatticeND",
    # Random
    "ErdosRenyi",
    "kRegularGraph",
    "ConfigurationModel",
    "RandomGeometric",
    "LFRBenchmark",
    # Scale-Free
    "BarabasiAlbert",
    "ExtendedBarabasiAlbert",
    "DualBarabasiAlbert",
    "HolmeKim",
    # Small-World
    "WattsStrogatz",
    # Community
    "StochasticBlockModel",
    # Complete
    "FullyConnected",
    # Bipartite
    "BipartiteGraph",
    "BipartiteFromDegreeSequence",
    # Neural
    "HofieldNN",
    "SCSGeneralizedNN",
    # Multispectral
    "MultiplicativeCascade",
    "VicsekGraph",
    "HierarchicalModular",
    # Fractal
    "SierpinskiGraph",
    "DGMgraph",
    # Graph-of-Graphs / Dirac
    "GraphOfGraphs",
    "DiracCombGraph",
    "DiracBrushGraph",
    # Temporal
    "TemporalGraph",
]
