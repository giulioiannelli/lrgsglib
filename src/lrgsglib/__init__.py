import sys
import warnings

from .shared import *
from .core import *
from .plotlib import *
from .config.const import DEFAULT_RECURSION_LIMIT, PATHDATA
from lrgsglib.config.funcs import move_to_rootf


# Explicit public API: graph types
from .graphs.nx import (
    Lattice2D,
    Lattice3D,
    SignedGraph,
    MultiplicativeCascadeGraph,
    SCSGeneralizedNN,
    load_or_compute_Lattice2D,
    load_or_compute_Lattice3D,
)

# Explicit public API: dynamics models
from .statsys import (
    IsingDynamics,
    ContactProcess,
    ContactProcessEI,
    ContactProcessSIR,
)

# Explicit public API: commonly used utilities
from .utils.basic import (
    flip_to_positive_majority,
    join_non_empty,
)
from .utils.tools import Chronometer

sys.setrecursionlimit(DEFAULT_RECURSION_LIMIT)
warnings.simplefilter(action="ignore", category=FutureWarning)

# Note: Logging is opt-in. Use lrgsglib.loglib.setup_custom_logger() if needed.
