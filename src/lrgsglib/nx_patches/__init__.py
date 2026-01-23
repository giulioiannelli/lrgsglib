# Base class
from .SignedGraph.SignedGraph import SignedGraph

# Legacy imports (backward compatibility)
from .ErdosRenyi import ErdosRenyi
from .FullyConnected import FullyConnected
from .HofieldNN import HofieldNN
from .SCSGeneralizedNN import SCSGeneralizedNN
from .DGMgraph import DGMgraph
from .Lattice2D import Lattice2D, load_or_compute_Lattice2D
from .Lattice3D import Lattice3D, load_or_compute_Lattice3D
from .MultispectralGraph import MultispectralGraph
from .MultiplicativeCascade import MultiplicativeCascadeGraph
from .HierarchicalModular import HierarchicalModularNetwork
from .Vicsek import VicsekGraph
from .DiracLattice import DiracLattice
from .WattStrogatz import WattStrogatz

# New module structure (unified hierarchy)
from .lattice import LatticeND
from .random import (
    RandomGraph, BarabasiAlbert, StochasticBlockModel,
    kRegularGraph, ConfigurationModel, RandomGeometric,
    LFRBenchmark, ExtendedBarabasiAlbert, DualBarabasiAlbert, HolmeKim,
)
from .complete import CompleteGraph
from .fractal import FractalGraph, SierpinskiGraph
from .bipartite import BipartiteGraph, BipartiteFromDegreeSequence
from .temporal import TemporalGraph, TemporalSignedGraph
from .datasets import RealDatasetLoader, BUILTIN_DATASETS

# Utility functions
from .funcs import *
