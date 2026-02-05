"""Random graph family module.

DEPRECATED: This module re-exports from lrgsglib.graphs.nx.
Import directly from lrgsglib.graphs.nx instead.
"""

import warnings

warnings.warn(
    "Importing from lrgsglib.nx_patches.random is deprecated. "
    "Use lrgsglib.graphs.nx instead.",
    DeprecationWarning,
    stacklevel=2
)

from ...graphs.nx.random import (
    RandomGraphNX as RandomGraph,
    RandomNwContainerBase,
    ErdosRenyiNX as ErdosRenyi,
    WattsStrogatzNX as WattsStrogatz,
    BarabasiAlbertNX as BarabasiAlbert,
    StochasticBlockModelNX as StochasticBlockModel,
    kRegularGraphNX as kRegularGraph,
    ConfigurationModelNX as ConfigurationModel,
    RandomGeometricNX as RandomGeometric,
    LFRBenchmarkNX as LFRBenchmark,
    ExtendedBarabasiAlbertNX as ExtendedBarabasiAlbert,
    DualBarabasiAlbertNX as DualBarabasiAlbert,
    HolmeKimNX as HolmeKim,
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
