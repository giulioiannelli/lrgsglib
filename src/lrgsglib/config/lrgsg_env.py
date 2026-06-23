"""Canonical filesystem locations for lrgsglib.

Each path resolves **env-first, else derived**:

* If the corresponding ``LRGSG_*`` environment variable is set (the dev
  activation script exports the full set), that value is used verbatim — so on a
  configured machine the values are exactly what they have always been.
* Otherwise the path is **derived from this file's location**, making the
  package relocatable (a fresh checkout or an installed wheel works with no
  environment setup). Code paths follow the *package*; output paths
  (``LRGSG_DATA``/``LRGSG_IPYNB``/``LRGSG_LOG``) follow the outer-repo root
  ``LRGSG_LLIB``.
"""

import os
from pathlib import Path


def _env(name: str, default: Path) -> Path:
    """Resolve ``LRGSG_<name>`` from the environment, else ``default``."""
    raw = os.environ.get(name)
    return Path(raw).resolve() if raw else default


# --- Anchors derived from this file ------------------------------------------
# .../<outer>/lrgsglib/src/lrgsglib/config/lrgsg_env.py
_PKG_DIR: Path = Path(__file__).resolve().parents[1]   # .../src/lrgsglib  (the package)
_SRC: Path = _PKG_DIR.parent                           # .../src
_SUBMODULE: Path = _SRC.parent                         # .../lrgsglib      (submodule root)
_OUTER: Path = _SUBMODULE.parent                       # .../<outer repo>  (contains data/)

_STATSYS: Path = _PKG_DIR / "statsys"
_GRAPHS: Path = _PKG_DIR / "graphs"
_GT: Path = _GRAPHS / "gt"


def _model_ccore(model: str) -> Path:
    """``ccore`` directory of a statsys model package (derivation default)."""
    return _STATSYS / model / "ccore"


# --- Outer-repo root + output/data dirs --------------------------------------
LRGSG_LLIB: Path = _env("LRGSG_LLIB", _OUTER)
LRGSG_DATA: Path = _env("LRGSG_DATA", LRGSG_LLIB / "data")
LRGSG_IPYNB: Path = _env("LRGSG_IPYNB", LRGSG_LLIB / "ipynb")
LRGSG_LOG: Path = _env("LRGSG_LOG", LRGSG_LLIB / ".log")

# --- Submodule-level code paths ----------------------------------------------
LRGSG_BUILD: Path = _env("LRGSG_BUILD", _SUBMODULE / "build")
LRGSG_SRC: Path = _env("LRGSG_SRC", _SRC)
LRGSG_TEST: Path = _env("LRGSG_TEST", _SUBMODULE / "test")
LRGSG_TOOLS: Path = _env("LRGSG_TOOLS", _SUBMODULE / "tools")
LRGSG_TOOLS_SCRPT: Path = _env("LRGSG_TOOLS_SCRPT", LRGSG_TOOLS / "bash")
LRGSG_TOOLS_PY: Path = _env("LRGSG_TOOLS_PY", LRGSG_TOOLS / "py")

# --- Package-level code paths ------------------------------------------------
LRGSG_LIB: Path = _env("LRGSG_LIB", _PKG_DIR)
LRGSG_LIB_STOCPROC: Path = _env("LRGSG_LIB_STOCPROC", _PKG_DIR / "stocproc")
LRGSG_LIB_STATSYS: Path = _env("LRGSG_LIB_STATSYS", _STATSYS)
LRGSG_STATSYS_CCORE: Path = _env("LRGSG_STATSYS_CCORE", _STATSYS / "_ccore")
LRGSG_CCORE_SFMT: Path = _env("LRGSG_CCORE_SFMT", LRGSG_STATSYS_CCORE / "SFMT")

# --- Per-model ccore trees ---------------------------------------------------
LRGSG_STATSYS_ISING: Path = _env("LRGSG_STATSYS_ISING", _model_ccore("IsingDynamics"))
LRGSG_STATSYS_ISING_BIN: Path = _env("LRGSG_STATSYS_ISING_BIN", LRGSG_STATSYS_ISING / "bin")
LRGSG_STATSYS_CP: Path = _env("LRGSG_STATSYS_CP", _model_ccore("ContactProcess"))
LRGSG_STATSYS_CP_BIN: Path = _env("LRGSG_STATSYS_CP_BIN", LRGSG_STATSYS_CP / "bin")
LRGSG_STATSYS_VM: Path = _env("LRGSG_STATSYS_VM", _model_ccore("VoterModel"))
LRGSG_STATSYS_VM_BIN: Path = _env("LRGSG_STATSYS_VM_BIN", LRGSG_STATSYS_VM / "bin")
LRGSG_STATSYS_SRW: Path = _env("LRGSG_STATSYS_SRW", _model_ccore("SignedRW"))
LRGSG_STATSYS_KUR: Path = _env("LRGSG_STATSYS_KUR", _model_ccore("KuramotoModel"))
LRGSG_STATSYS_KUR_BIN: Path = _env("LRGSG_STATSYS_KUR_BIN", LRGSG_STATSYS_KUR / "bin")
LRGSG_STATSYS_RD: Path = _env("LRGSG_STATSYS_RD", _model_ccore("ReactionDiffusionModel"))
LRGSG_STATSYS_RD_BIN: Path = _env("LRGSG_STATSYS_RD_BIN", LRGSG_STATSYS_RD / "bin")
LRGSG_STATSYS_CODE: Path = _env("LRGSG_STATSYS_CODE", _model_ccore("CoupledODEModel"))
LRGSG_STATSYS_CODE_BIN: Path = _env("LRGSG_STATSYS_CODE_BIN", LRGSG_STATSYS_CODE / "bin")
LRGSG_STATSYS_POTTS: Path = _env("LRGSG_STATSYS_POTTS", _model_ccore("PottsModel"))
LRGSG_STATSYS_POTTS_BIN: Path = _env("LRGSG_STATSYS_POTTS_BIN", LRGSG_STATSYS_POTTS / "bin")
LRGSG_STATSYS_XY: Path = _env("LRGSG_STATSYS_XY", _model_ccore("XYModel"))
LRGSG_STATSYS_XY_BIN: Path = _env("LRGSG_STATSYS_XY_BIN", LRGSG_STATSYS_XY / "bin")
LRGSG_STATSYS_HBERG: Path = _env("LRGSG_STATSYS_HBERG", _model_ccore("HeisenbergModel"))
LRGSG_STATSYS_HBERG_BIN: Path = _env("LRGSG_STATSYS_HBERG_BIN", LRGSG_STATSYS_HBERG / "bin")
LRGSG_STATSYS_MSPEC: Path = _env("LRGSG_STATSYS_MSPEC", _model_ccore("MultiSpeciesModel"))
LRGSG_STATSYS_MSPEC_BIN: Path = _env("LRGSG_STATSYS_MSPEC_BIN", LRGSG_STATSYS_MSPEC / "bin")

# --- Ising (RBIM) C++ subtrees & SignedRW lattices ---------------------------
LRGSG_RBIM_BASE: Path = _env("LRGSG_RBIM_BASE", LRGSG_STATSYS_ISING / "base")
LRGSG_RBIM_STORE: Path = _env("LRGSG_RBIM_STORE", LRGSG_STATSYS_ISING / "storer")
LRGSG_RBIM_NATIVE: Path = _env("LRGSG_RBIM_NATIVE", LRGSG_STATSYS_ISING / "native")
LRGSG_SRW_LATT: Path = _env("LRGSG_SRW_LATT", LRGSG_STATSYS_SRW / "Lattices")

# --- Graph engine code paths -------------------------------------------------
LRGSG_LIB_GRAPHS: Path = _env("LRGSG_LIB_GRAPHS", _GRAPHS)
LRGSG_GRAPHS_GT: Path = _env("LRGSG_GRAPHS_GT", _GT)
LRGSG_GRAPHS_GT_L2D: Path = _env("LRGSG_GRAPHS_GT_L2D", _GT / "Lattice2DGT")
LRGSG_GRAPHS_GT_L2D_CPP: Path = _env("LRGSG_GRAPHS_GT_L2D_CPP", LRGSG_GRAPHS_GT_L2D / "cpp")
LRGSG_GRAPHS_GT_BA: Path = _env("LRGSG_GRAPHS_GT_BA", _GT / "BarabasiAlbertGT")
LRGSG_GRAPHS_GT_BA_CPP: Path = _env("LRGSG_GRAPHS_GT_BA_CPP", LRGSG_GRAPHS_GT_BA / "cpp")
LRGSG_GRAPHS_GT_EBA: Path = _env("LRGSG_GRAPHS_GT_EBA", _GT / "ExtendedBarabasiAlbertGT")
LRGSG_GRAPHS_GT_EBA_CPP: Path = _env("LRGSG_GRAPHS_GT_EBA_CPP", LRGSG_GRAPHS_GT_EBA / "cpp")
LRGSG_GRAPHS_GT_DBA: Path = _env("LRGSG_GRAPHS_GT_DBA", _GT / "DualBarabasiAlbertGT")
LRGSG_GRAPHS_GT_DBA_CPP: Path = _env("LRGSG_GRAPHS_GT_DBA_CPP", LRGSG_GRAPHS_GT_DBA / "cpp")
LRGSG_GRAPHS_GT_FC: Path = _env("LRGSG_GRAPHS_GT_FC", _GT / "FullyConnectedGT")
LRGSG_GRAPHS_GT_FC_CPP: Path = _env("LRGSG_GRAPHS_GT_FC_CPP", LRGSG_GRAPHS_GT_FC / "cpp")
LRGSG_GRAPHS_GT_HK: Path = _env("LRGSG_GRAPHS_GT_HK", _GT / "HolmeKimGT")
LRGSG_GRAPHS_GT_HK_CPP: Path = _env("LRGSG_GRAPHS_GT_HK_CPP", LRGSG_GRAPHS_GT_HK / "cpp")
LRGSG_GRAPHS_GT_MC: Path = _env("LRGSG_GRAPHS_GT_MC", _GT / "MultiplicativeCascadeGT")
LRGSG_GRAPHS_GT_MC_CPP: Path = _env("LRGSG_GRAPHS_GT_MC_CPP", LRGSG_GRAPHS_GT_MC / "cpp")
