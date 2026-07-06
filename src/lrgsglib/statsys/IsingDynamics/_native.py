"""Loader + protocol constants for the ``_ising_native`` pybind11 module.

One shared home for the pieces every pb-wired scheme class needs: the
import-with-a-clear-error loader (the legacy god-object carries its own
private copy; new code imports from here) and the mapping from the Layer-2
``order`` axis to the update-mode strings the compiled kernel understands
(``LRGSG_rbim.c::glauberMetropolis_Nstep``: ``'sequential'`` = fixed 0..N-1
order, ``'asynchronous'`` = uniform random sites with replacement). The
python axis value ``'permutation'`` has no native counterpart — schemes
reject it at ``supports()`` time.
"""

from __future__ import annotations

from typing import Any

__all__ = ["load_ising_native", "PB_ORDER_MAP"]

#: Layer-2 ``order`` axis value -> native kernel update-mode string.
PB_ORDER_MAP: dict[str, str] = {
    "random": "asynchronous",
    "typewriter": "sequential",
}


def load_ising_native() -> Any:
    """Import the compiled ``_ising_native`` module or fail loudly."""
    try:
        from .ccore import _ising_native  # type: ignore[import-untyped]
    except ImportError as exc:
        raise RuntimeError(
            "Pybind11 Ising backend not available. Build the extension with "
            "`make` in IsingDynamics/ccore/native/ (or `pip install -e .`)."
        ) from exc
    return _ising_native
