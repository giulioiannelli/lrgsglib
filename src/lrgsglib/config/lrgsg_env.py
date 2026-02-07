"""Runtime path detection for lrgsglib.

All paths are resolved relative to the package location at import time.
Environment variables override defaults for DATA, IPYNB, LOG, and CCORE_BIN.

This replaces the former Makefile-generated file with hardcoded absolute paths.
"""

from pathlib import Path
import os

# ---------------------------------------------------------------------------
# Anchor points
# ---------------------------------------------------------------------------
# Package root: lrgsglib/src/lrgsglib/
_PKG_ROOT = Path(__file__).resolve().parent.parent

# Library root: lrgsglib/ (the submodule / repo root)
_LIB_ROOT = _PKG_ROOT.parent.parent

# Detect outer project root (if used as submodule in lrgsglib-ipynb)
# Heuristic: look for ipynb/ directory in parent of library root
_PROJECT_ROOT = (
    _LIB_ROOT.parent
    if (_LIB_ROOT.parent / "ipynb").exists()
    else _LIB_ROOT
)


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------
def _path(env_var: str, default: Path) -> Path:
    """Return Path from environment variable, falling back to default."""
    val = os.environ.get(env_var)
    return Path(val) if val else default


# ---------------------------------------------------------------------------
# Library structure paths (relative to package/library root)
# ---------------------------------------------------------------------------
LRGSG_LIB: Path = _PKG_ROOT
LRGSG_BUILD: Path = _LIB_ROOT / "build"
LRGSG_SRC: Path = _LIB_ROOT / "src"
LRGSG_TEST: Path = _LIB_ROOT / "test"
LRGSG_TOOLS: Path = _LIB_ROOT / "tools"
LRGSG_TOOLS_SCRPT: Path = _LIB_ROOT / "tools" / "bash"
LRGSG_TOOLS_PY: Path = _LIB_ROOT / "tools" / "py"

# ---------------------------------------------------------------------------
# Package subpackage paths
# ---------------------------------------------------------------------------
LRGSG_LIB_CCORE: Path = _PKG_ROOT / "Ccore"
LRGSG_LIB_GT_PATCHES: Path = _PKG_ROOT / "gt_patches"
LRGSG_LIB_NX_PATCHES: Path = _PKG_ROOT / "nx_patches"
LRGSG_LIB_STOCPROC: Path = _PKG_ROOT / "stocproc"

# ---------------------------------------------------------------------------
# Ccore paths
# ---------------------------------------------------------------------------
LRGSG_CCORE_BIN: Path = _path("LRGSG_CCORE_BIN", _PKG_ROOT / "Ccore" / "bin")
LRGSG_CCORE_SFMT: Path = _PKG_ROOT / "Ccore" / "SFMT"
LRGSG_CCORE_STATSYS: Path = _PKG_ROOT / "Ccore" / "statsys"
LRGSG_GT_PATCHES_CPP: Path = _PKG_ROOT / "gt_patches" / "cpp"

# ---------------------------------------------------------------------------
# Dynamics subsystem paths
# ---------------------------------------------------------------------------
LRGSG_STATSYS_RBIM: Path = _PKG_ROOT / "Ccore" / "statsys" / "RBIsingM"
LRGSG_STATSYS_SRW: Path = _PKG_ROOT / "Ccore" / "statsys" / "signedRw"
LRGSG_STATSYS_VM: Path = _PKG_ROOT / "Ccore" / "statsys" / "voterM"
LRGSG_STATSYS_CP: Path = _PKG_ROOT / "Ccore" / "statsys" / "contactP"
LRGSG_RBIM_BASE: Path = _PKG_ROOT / "Ccore" / "statsys" / "RBIsingM" / "base"
LRGSG_RBIM_SIMC: Path = (
    _PKG_ROOT / "Ccore" / "statsys" / "RBIsingM" / "simulatorC"
)
LRGSG_RBIM_STORE: Path = (
    _PKG_ROOT / "Ccore" / "statsys" / "RBIsingM" / "storer"
)
LRGSG_SRW_LATT: Path = (
    _PKG_ROOT / "Ccore" / "statsys" / "signedRw" / "Lattices"
)

# ---------------------------------------------------------------------------
# User-facing / project-level paths (overridable via environment)
# ---------------------------------------------------------------------------
LRGSG_DATA: Path = _path("LRGSG_DATA", _PROJECT_ROOT / "data")
LRGSG_IPYNB: Path = _path("LRGSG_IPYNB", _PROJECT_ROOT / "ipynb")
LRGSG_LOG: Path = _path("LRGSG_LOG", _PROJECT_ROOT / ".log")

# ---------------------------------------------------------------------------
# Outer project root
# ---------------------------------------------------------------------------
LRGSG_LLIB: Path = _PROJECT_ROOT
