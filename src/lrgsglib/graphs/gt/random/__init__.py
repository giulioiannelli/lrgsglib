"""
graph-tool random graph implementations.

This module provides random graph classes with graph-tool backend.
Classes without native GT implementations fall back to NetworkX.
"""

from .ErdosRenyiGT import ErdosRenyiGT
from .BarabasiAlbertGT import BarabasiAlbertGT
from .WattsStrogatzGT import WattsStrogatzGT
from .StochasticBlockModelGT import StochasticBlockModelGT
from .RandomGraphGT import RandomGraphGT
from .kRegularGraphGT import kRegularGraphGT
from .ConfigurationModelGT import ConfigurationModelGT
from .RandomGeometricGT import RandomGeometricGT
from .LFRBenchmarkGT import LFRBenchmarkGT
from .ExtendedBarabasiAlbertGT import ExtendedBarabasiAlbertGT
from .DualBarabasiAlbertGT import DualBarabasiAlbertGT
from .HolmeKimGT import HolmeKimGT

__all__ = [
    # Native GT implementations
    "ErdosRenyiGT",
    "BarabasiAlbertGT",
    "WattsStrogatzGT",
    "StochasticBlockModelGT",
    # Fallback to NX
    "RandomGraphGT",
    "kRegularGraphGT",
    "ConfigurationModelGT",
    "RandomGeometricGT",
    "LFRBenchmarkGT",
    "ExtendedBarabasiAlbertGT",
    "DualBarabasiAlbertGT",
    "HolmeKimGT",
]
