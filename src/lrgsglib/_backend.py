"""Package-level home for the numerical array-backend manager.

Re-exports the proven numpy/scipy/cupy dispatcher implemented in
``graphs/_shared/_backend.py`` so that *both* dynamics solvers (``statsys``) and
spectral utilities (``utils.lrg``) can share one GPU/CPU backend selector
without reaching into the graphs package internals. Importing from here keeps
the dependency direction clean (consumers depend on ``lrgsglib._backend``, not
on a sibling subpackage's private module).
"""

from .graphs._shared._backend import ArrayBackend, Backend, BackendManager

__all__ = ["ArrayBackend", "Backend", "BackendManager"]
