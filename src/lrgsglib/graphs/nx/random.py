"""Random graph facade module.

Re-exports random graph classes from their respective packages for cleaner imports.
"""

from .BarabasiAlbertNX import (  # Backward compatibility alias
    BarabasiAlbert,
    BarabasiAlbertNX,
)
from .ConfigurationModelNX import (  # Backward compatibility alias
    ConfigurationModel,
    ConfigurationModelNX,
)
from .DualBarabasiAlbertNX import (  # Backward compatibility alias
    DualBarabasiAlbert,
    DualBarabasiAlbertNX,
)
from .ErdosRenyiNX import (  # Backward compatibility alias
    ErdosRenyi,
    ErdosRenyiNX,
)
from .ExtendedBarabasiAlbertNX import (  # Backward compatibility alias
    ExtendedBarabasiAlbert,
    ExtendedBarabasiAlbertNX,
)
from .HolmeKimNX import (  # Backward compatibility alias
    HolmeKim,
    HolmeKimNX,
)
from .kRegularGraphNX import (  # Backward compatibility alias
    kRegularGraph,
    kRegularGraphNX,
)
from .LFRBenchmarkNX import (  # Backward compatibility alias
    LFRBenchmark,
    LFRBenchmarkNX,
)
from .RandomGeometricNX import (  # Backward compatibility alias
    RandomGeometric,
    RandomGeometricNX,
)
from .RandomGraphNX import (  # Backward compatibility alias
    RandomGraph,
    RandomGraphNX,
    RandomNwContainerBase,
)
from .StochasticBlockModelNX import (  # Backward compatibility alias
    StochasticBlockModel,
    StochasticBlockModelNX,
)
from .WattsStrogatzNX import (  # Backward compatibility alias
    WattsStrogatz,
    WattsStrogatzNX,
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
