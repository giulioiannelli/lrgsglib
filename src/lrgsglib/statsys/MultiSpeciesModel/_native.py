"""Loader for the in-process multi-species pybind11 kernel
(``_multispec_native``).

The compiled module wraps the SAME C kernels the legacy subprocess backend
uses (``LRGSG_multispec.c``), re-seeded deterministically per call, so pb
runs are seed-reproducible. The kernel implements ONE update: Metropolis
acceptance, uniform-random site visits WITH replacement, a uniform random
(component, other-label) proposal, every ``ΔE <= 0`` accepted — and it
assumes the IDENTITY interaction matrix and a UNIFORM number of states per
species. The per-state extensive energy is recorded AFTER each sweep (energy
only — no order parameter). Anything else is a python-backend capability
(``MultiSpeciesMetropolis._pb_check_supported`` enforces the boundary).
"""

from __future__ import annotations

from types import ModuleType

__all__ = ["load_multispec_native"]


def load_multispec_native() -> ModuleType:
    """Import and return ``.ccore._multispec_native``, or raise with a hint."""
    try:
        from .ccore import _multispec_native
    except ImportError as exc:
        raise RuntimeError(
            "The native multi-species kernel (_multispec_native) is not "
            "built. Build it with `make` in "
            "src/lrgsglib/statsys/MultiSpeciesModel/ccore/native/ (needs the "
            "_ccore/SFMT submodule initialised), or use runlang='py'."
        ) from exc
    return _multispec_native
