"""
Random graph family module.

This module provides a clean class hierarchy for random graphs:
- RandomGraph: Abstract base class for random graphs
- ErdosRenyi: Erdos-Renyi random graph (giant component)
- WattsStrogatz: Small-world graph
- BarabasiAlbert: Scale-free graph (preferential attachment)
- StochasticBlockModel: Community-structured random graph
- kRegularGraph: Random k-regular graph
- ConfigurationModel: Graph with prescribed degree sequence
- RandomGeometric: Spatial random graph
- LFRBenchmark: Community detection benchmark graph
- ExtendedBarabasiAlbert: BA with initial attractiveness
- DualBarabasiAlbert: BA with two attachment modes
- HolmeKim: BA with clustering (triad formation)

Example
-------
>>> from lrgsglib.nx_patches.random import ErdosRenyi, BarabasiAlbert
>>> er = ErdosRenyi(n=100, p=0.05, pflip=0.1)
>>> ba = BarabasiAlbert(n=100, m=3, pflip=0.1)
"""

from .random_graph import RandomGraph, RandomNwContainerBase
from .erdos_renyi import ErdosRenyi
from .watts_strogatz import WattsStrogatz
from .barabasi_albert import BarabasiAlbert
from .stochastic_block import StochasticBlockModel
from .k_regular import kRegularGraph
from .configuration_model import ConfigurationModel
from .random_geometric import RandomGeometric
from .lfr_benchmark import LFRBenchmark
from .preferential_attachment import (
    ExtendedBarabasiAlbert,
    DualBarabasiAlbert,
    HolmeKim,
)

__all__ = [
    "RandomGraph",
    "RandomNwContainerBase",
    "ErdosRenyi",
    "WattsStrogatz",
    "BarabasiAlbert",
    "StochasticBlockModel",
    "kRegularGraph",
    "ConfigurationModel",
    "RandomGeometric",
    "LFRBenchmark",
    "ExtendedBarabasiAlbert",
    "DualBarabasiAlbert",
    "HolmeKim",
]
