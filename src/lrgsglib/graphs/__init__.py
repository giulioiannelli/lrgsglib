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
Lattice2D : 2D lattice (square, triangular, hexagonal)
    - Parameters: side1, side2, geo, pflip, periodic, seed, engine

Lattice3D : 3D lattice (simple cubic, BCC, FCC)
    - Parameters: dim, geo, pflip, periodic, seed, engine

ErdosRenyi : Erdos-Renyi random graph
    - Parameters: n, p, pflip, extract_giant_component, seed, engine

BarabasiAlbert : Barabasi-Albert scale-free graph
    - Parameters: n, m, pflip, seed, engine

WattsStrogatz : Watts-Strogatz small-world graph
    - Parameters: n, k, p, pflip, seed, engine

StochasticBlockModel : Stochastic Block Model with communities
    - Parameters: sizes, p_matrix, pflip, extract_giant_component, seed, engine

Engine Configuration
--------------------
Engines can be configured via:
1. Explicit parameter: `Lattice2D(..., engine='gt')`
2. Global setting: `set_default_engine('gt')`
3. Environment variable: `LRGSG_GRAPH_ENGINE=gt`

Priority: explicit parameter > environment variable > global setting

Supported Engines
----------------
- 'nx' : NetworkX (default, always available)
- 'gt' : graph-tool (high performance, requires installation)
- 'ig' : igraph (future support)

Module Structure
----------------
The module is organized as follows:

graphs/
├── __init__.py              # This file (unified API)
├── protocols.py             # Protocol definitions
├── _engine.py               # Engine selection
├── Lattice2D.py             # Facade
├── Lattice3D.py             # Facade
├── ErdosRenyi.py            # Facade
├── BarabasiAlbert.py        # Facade
├── WattsStrogatz.py         # Facade
├── StochasticBlockModel.py  # Facade
├── nx/                      # NetworkX implementations
│   ├── base/SignedGraphNX.py
│   ├── lattice/{Lattice2DNX,Lattice3DNX}.py
│   └── random/{ErdosRenyiNX,...}.py
└── gt/                      # graph-tool implementations
    ├── base/SignedGraphGT.py
    ├── lattice/{Lattice2DGT,Lattice3DGT}.py
    └── random/{ErdosRenyiGT,...}.py

Direct access to engine-specific implementations:
>>> from lrgsglib.graphs.nx.lattice import Lattice2DNX
>>> from lrgsglib.graphs.gt.lattice import Lattice2DGT
"""

from .protocols import (
    LatticeGraphProtocol,
    SignedGraphProtocol,
    SpectralGraphProtocol,
    is_lattice_graph,
    is_signed_graph,
)
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

# Import graph type facades from root level
from .Lattice2D import Lattice2D
from .Lattice3D import Lattice3D
from .ErdosRenyi import ErdosRenyi
from .BarabasiAlbert import BarabasiAlbert
from .WattsStrogatz import WattsStrogatz
from .StochasticBlockModel import StochasticBlockModel
from .GraphOfGraphs import GraphOfGraphs

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
    # Graph types - Lattice
    "Lattice2D",
    "Lattice3D",
    # Graph types - Random
    "ErdosRenyi",
    "BarabasiAlbert",
    "WattsStrogatz",
    "StochasticBlockModel",
    # Graph types - Hierarchical
    "GraphOfGraphs",
]
