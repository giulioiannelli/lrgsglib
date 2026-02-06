"""Random graph facade module.

Re-exports random graph classes from their respective packages for cleaner imports.
"""

from .RandomGraphNX import (
    RandomGraphNX,
    RandomNwContainerBase,
    # Backward compatibility alias
    RandomGraph,
)
from .ErdosRenyiNX import (
    ErdosRenyiNX,
    # Backward compatibility alias
    ErdosRenyi,
)
from .BarabasiAlbertNX import (
    BarabasiAlbertNX,
    # Backward compatibility alias
    BarabasiAlbert,
)
from .WattsStrogatzNX import (
    WattsStrogatzNX,
    # Backward compatibility alias
    WattsStrogatz,
)
from .StochasticBlockModelNX import (
    StochasticBlockModelNX,
    # Backward compatibility alias
    StochasticBlockModel,
)
from .kRegularGraphNX import (
    kRegularGraphNX,
    # Backward compatibility alias
    kRegularGraph,
)
from .ConfigurationModelNX import (
    ConfigurationModelNX,
    # Backward compatibility alias
    ConfigurationModel,
)
from .RandomGeometricNX import (
    RandomGeometricNX,
    # Backward compatibility alias
    RandomGeometric,
)
from .LFRBenchmarkNX import (
    LFRBenchmarkNX,
    # Backward compatibility alias
    LFRBenchmark,
)
from .ExtendedBarabasiAlbertNX import (
    ExtendedBarabasiAlbertNX,
    # Backward compatibility alias
    ExtendedBarabasiAlbert,
)
from .DualBarabasiAlbertNX import (
    DualBarabasiAlbertNX,
    # Backward compatibility alias
    DualBarabasiAlbert,
)
from .HolmeKimNX import (
    HolmeKimNX,
    # Backward compatibility alias
    HolmeKim,
)

__all__ = [
    # Base
    "RandomGraphNX",
    "RandomGraph",
    "RandomNwContainerBase",
    # Erdos-Renyi
    "ErdosRenyiNX",
    "ErdosRenyi",
    # Barabasi-Albert
    "BarabasiAlbertNX",
    "BarabasiAlbert",
    # Watts-Strogatz
    "WattsStrogatzNX",
    "WattsStrogatz",
    # Stochastic Block Model
    "StochasticBlockModelNX",
    "StochasticBlockModel",
    # k-Regular
    "kRegularGraphNX",
    "kRegularGraph",
    # Configuration Model
    "ConfigurationModelNX",
    "ConfigurationModel",
    # Random Geometric
    "RandomGeometricNX",
    "RandomGeometric",
    # LFR Benchmark
    "LFRBenchmarkNX",
    "LFRBenchmark",
    # Preferential Attachment Variants
    "ExtendedBarabasiAlbertNX",
    "ExtendedBarabasiAlbert",
    "DualBarabasiAlbertNX",
    "DualBarabasiAlbert",
    "HolmeKimNX",
    "HolmeKim",
]
