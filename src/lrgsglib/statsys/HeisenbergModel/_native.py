"""Loader for the in-process Heisenberg pybind11 kernel (``_heisenberg_native``).

The compiled module wraps the SAME C kernels the legacy subprocess backend
uses (``LRGSG_heisenberg.c``), re-seeded deterministically per call, so pb
runs are seed-reproducible. The kernel implements ONE update: Metropolis
acceptance, uniform-random site visits WITH replacement, proposal
``normalize(n + delta·u)`` with ``u`` uniform in [−1, 1]³, every ``ΔE <= 0``
accepted; observables (extensive pairwise energy + magnetisation
``|Σ n_i|/N``) are recorded AFTER each sweep. Anything else is a
python-backend capability (``HeisenbergMetropolis._pb_check_supported``
enforces the boundary).
"""

from __future__ import annotations

from types import ModuleType

__all__ = ["load_heisenberg_native"]


def load_heisenberg_native() -> ModuleType:
    """Import and return ``.ccore._heisenberg_native``, or raise with a hint."""
    try:
        from .ccore import _heisenberg_native
    except ImportError as exc:
        raise RuntimeError(
            "The native Heisenberg kernel (_heisenberg_native) is not built. "
            "Build it with `make` in "
            "src/lrgsglib/statsys/HeisenbergModel/ccore/native/ (needs the "
            "_ccore/SFMT submodule initialised), or use runlang='py'."
        ) from exc
    return _heisenberg_native
