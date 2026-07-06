"""Loader for the in-process XY pybind11 kernel (``_xy_native``).

The compiled module wraps the SAME C kernels the legacy subprocess backend
uses (``LRGSG_xy.c``), re-seeded deterministically per call, so pb runs are
seed-reproducible. The kernel implements ONE update: Metropolis acceptance,
uniform-random site visits WITH replacement, proposal uniform in
``[θ − delta, θ + delta]`` (mod 2π), every ``ΔE <= 0`` accepted; observables
(extensive pairwise energy + magnetisation ``|Σ e^{iθ}|/N``) are recorded
AFTER each sweep. Anything else is a python-backend capability
(``XYMetropolis._pb_check_supported`` enforces the boundary).
"""

from __future__ import annotations

from types import ModuleType

__all__ = ["load_xy_native"]


def load_xy_native() -> ModuleType:
    """Import and return ``.ccore._xy_native``, or raise with a build hint."""
    try:
        from .ccore import _xy_native
    except ImportError as exc:
        raise RuntimeError(
            "The native XY kernel (_xy_native) is not built. Build it with "
            "`make` in src/lrgsglib/statsys/XYModel/ccore/native/ (needs the "
            "_ccore/SFMT submodule initialised), or use runlang='py'."
        ) from exc
    return _xy_native
