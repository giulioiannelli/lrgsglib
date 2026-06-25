"""Default constants for the SignedRW dynamics family.

The package ships two independent dynamics models, each dispatching through the
shared solver registry under its own name:

* the single-walker family (``SignedWalker`` / ``Absorbing``/``Killing``/
  ``StickyWalker``): a ``py`` reference kernel + a ``pb_absorb`` pybind11 kernel;
* the spin-copy-with-sign dynamics (``SignedSpinCopy``, legacy alias
  ``SignedRW``): a ``py`` loop + a ``C`` subprocess.

The registry names mirror each class's ``dyn_UVclass`` identifier.
"""

from __future__ import annotations

__all__ = [
    "SRW_WALKER_SOLVER_NAME",
    "SRW_COPY_SOLVER_NAME",
]

# Solver-registry names (keys in ``statsys._solver_engine``).
SRW_WALKER_SOLVER_NAME: str = "signed_walker"
SRW_COPY_SOLVER_NAME: str = "signed_spin_copy"
