import sys
import warnings

# Package version. setuptools_scm writes ``_version.py`` at build time
# (single source of truth = annotated git tags ``vX.Y.Z``). Fall back to the
# installed distribution metadata, then to a sentinel, so ``__version__`` is
# always defined at runtime even in an uninstalled / source checkout.
try:
    from ._version import version as __version__
except Exception:  # pragma: no cover - depends on build/install state
    from importlib.metadata import version as _pkg_version

    try:
        __version__ = _pkg_version("lrgsglib")
    except Exception:  # pragma: no cover - not installed as a distribution
        __version__ = "0.0.0+unknown"

from lrgsglib.config.funcs import move_to_rootf

from .config.const import DEFAULT_RECURSION_LIMIT, PATHDATA
from .core import *

# Explicit public API: graph types
from .graphs.nx import (
    Lattice2D,
    Lattice3D,
    MultiplicativeCascadeGraph,
    SCSGeneralizedNN,
    SignedGraph,
    load_or_compute_Lattice2D,
    load_or_compute_Lattice3D,
)
from .plotlib import *
from .shared import (  # Core libraries users expect at top level; Standard-library modules used by scripts/serializers; Common types; Commonly used constants; Frequently used scipy.cluster / scipy.signal helpers
    Decimal,
    Fraction,
    Graph,
    Iterable,
    NDArray,
    Number,
    Path,
    cophenet,
    copy,
    cp,
    dendrogram,
    expm,
    fcluster,
    find_peaks,
    leaves_list,
    linkage,
    np,
    nx,
    os,
    pd,
    pdist,
    pk,
    plt,
    random,
    re,
    scipy,
    squareform,
    subprocess,
    sys,
    time,
    tqdm,
    warnings,
)

# Explicit public API: dynamics models
from .statsys import (
    ContactProcess,
    ContactProcessEI,
    ContactProcessSIR,
    IsingDynamics,
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
