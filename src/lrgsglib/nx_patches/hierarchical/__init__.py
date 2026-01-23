"""
Hierarchical/multispectral graph family module.

This module provides a clean interface for hierarchical and multispectral
graph types, which share a common abstract base class (MultispectralGraph).

Graph Types
-----------
- MultispectralGraph: Abstract base class
- MultiplicativeCascadeGraph: Scale-free hierarchical cascading structure
- VicsekGraph: Vicsek fractal graph
- HierarchicalModularNetwork: Modular hierarchical network
- DiracLatticeGraph: Dirac lattice base class
  - DiracCombGraph: Comb-shaped Dirac structure
  - DiracBrushGraph: Brush-shaped Dirac structure

Example
-------
>>> from lrgsglib.nx_patches.hierarchical import (
...     MultiplicativeCascadeGraph,
...     VicsekGraph,
...     HierarchicalModularNetwork,
... )
>>> mc = MultiplicativeCascadeGraph(p1=0.8, p2=0.6, iterations=7, pflip=0.1)
>>> mc.flip_random_fract_edges()
"""

# Re-export from existing modules (maintains backward compatibility)
from ..MultispectralGraph.MultispectralGraph import MultispectralGraph
from ..MultiplicativeCascade.MultiplicativeCascade import MultiplicativeCascadeGraph
from ..Vicsek.Vicsek import VicsekGraph
from ..HierarchicalModular.HierarchicalModular import HierarchicalModularNetwork
from ..DiracLattice.DiracLattice import DiracLatticeGraph
from ..DiracLattice.DiracComb import DiracCombGraph
from ..DiracLattice.DiracBrush import DiracBrushGraph

__all__ = [
    # Base class
    "MultispectralGraph",
    # Concrete implementations
    "MultiplicativeCascadeGraph",
    "VicsekGraph",
    "HierarchicalModularNetwork",
    # Dirac family
    "DiracLatticeGraph",
    "DiracCombGraph",
    "DiracBrushGraph",
]
