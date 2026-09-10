"""Session setup: root move, mpl style, paths, rng, verbosity, constants.

These helpers back the canonical lab/report header described in
``.agents/rules/lab-hygiene.md`` and ``.agents/rules/notebook-hygiene.md``::

    from lrgsglib.notebooks import *
    use_lab_style()
    paths = make_lab_paths("<CODE>")
    seed  = None
    rng   = resolved_rng(seed)

The ``nprint`` / ``nlog`` shims respect package-level flags so diagnostic
output stays opt-in (default OFF). The flags are read off the
``lrgsglib.notebooks`` package, not this module, so the documented
pattern keeps working::

    import lrgsglib.notebooks as _nb
    _nb.verbose_print_nb = True
"""
import sys as _sys
import logging as _logging
from pathlib import Path as _Path
import numpy as _np

from ..config.funcs import (
    move_to_rootf,
    peq_fstr,
    Teq_fstr,
    avgeq_parse,
    build_fname_or_pattern_direct,
    build_pT_fname,
)

# Most-used library defaults (symbolic, reusable).
from ..config.const import (
    PATHDATA,
    PATHPLOT,
    LRGSG_SRC,
    LRSG_ENTROPY_STEP,
    DEFAULT_ENTROPY_LEXPONENT,
    DEFAULT_ENTROPY_HEXPONENT,
    L2D_SIDE1,
    L2D_GEO_SQR,
    L2D_GEO_TRI,
    L2D_P_C_DICT,
    # CHL-SF04 protein-TMD experiment grid.
    PFLIP_SWEEP_GRID,
    K_FRAC_GALLERY,
    K_FRAC_Q3,
    N_SEEDS_PFLIP,
    CHL_SF04_MASTER_SEED,
    CHL_SF04_RES_MIN,
    CHL_SF04_RES_MAX,
    CHL_SF04_CORPUS_SIZE,
    CHL_SF04_N_FEATURES,
    CHL_SF04_LOG_MSE_FLOOR,
    Q3_EMERGENCE_THRESHOLD,
)

verbose_print_nb: bool = False
verbose_log_nb: bool = False


def _pkg_flag(name: str) -> bool:
    """Read a verbosity flag from the ``lrgsglib.notebooks`` package.

    Falls back to this module's own default when the package is not yet
    fully initialised (e.g. during import).
    """
    pkg = _sys.modules.get(__package__)
    return bool(getattr(pkg, name, globals()[name]))


def nprint(*args, **kwargs) -> None:
    """`print` shim: emits only when ``lrgsglib.notebooks.verbose_print_nb`` is True.

    Flip the flag in one cell to enable loud diagnostics within a lab
    or notebook session; leave it off everywhere else. Default OFF.
    """
    if _pkg_flag("verbose_print_nb"):
        print(*args, **kwargs)


def nlog(msg, level: int = _logging.INFO, *args, **kwargs) -> None:
    """`logging` shim: emits only when ``lrgsglib.notebooks.verbose_log_nb`` is True.

    Logs under the ``lrgsg.nb`` logger. Default OFF.
    """
    if _pkg_flag("verbose_log_nb"):
        _logging.getLogger("lrgsg.nb").log(level, msg, *args, **kwargs)


def resolved_rng(seed=None) -> "_np.random.Generator":
    """Return ``np.random.default_rng(seed)``. ``None`` → fresh entropy; int → reproducible.

    Use this exactly once per lab header so every stochastic call threads the
    same generator. Set ``seed`` to an int only when reproducing a specific
    result — hardcoding seeds in iterative exploration introduces systematic
    bias.
    """
    return _np.random.default_rng(seed)


def make_lab_paths(code_id: str) -> dict:
    """Canonical path dict for a project ``<code_id>``.

    Returns four sub-paths rooted at ``data/<code_id>/``:

        {"raw":     Path("data/<code_id>/raw"),
         "figures": Path("data/<code_id>/figures"),
         "cache":   Path("data/<code_id>/cache"),
         "cfg":     Path("data/<code_id>/cfg")}

    All parent directories are created on demand (``mkdir(parents=True,
    exist_ok=True)``). The returned dict is fresh on every call — callers
    may extend it with per-lab overrides.
    """
    root = _Path("data") / code_id
    out = {
        "raw":     root / "raw",
        "figures": root / "figures",
        "cache":   root / "cache",
        "cfg":     root / "cfg",
    }
    for p in out.values():
        p.mkdir(parents=True, exist_ok=True)
    return out


def use_lab_style() -> None:
    """Apply the interactive-environment mpl style sheet (``ipy/nb_plotsheet.mplstyle``).

    Safe to call multiple times; idempotent. Requires the current working
    directory to be the outer-repo root (``move_to_rootf`` handles this
    at package-import time).
    """
    import matplotlib.pyplot as _plt
    _plt.style.use("ipy/nb_plotsheet.mplstyle")
