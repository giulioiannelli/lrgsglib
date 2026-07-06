"""Loader for the in-process pybind11 Potts kernel (``_potts_native``).

One import point shared by the new scheme classes and the solver registry's
availability probe. The kernel (``ccore/LRGSG_potts.c`` via
``ccore/potts_native.cpp``) implements ONE update: Metropolis acceptance,
uniform random site visits WITH replacement, uniform random other-label
proposals, ties (ΔE = 0) always accepted — at every temperature, including
T = 0 (its acceptance test is ``dE <= 0`` first). Everything the kernel
cannot represent is a hard capability error in the schemes'
``_pb_check_supported`` (invariant #3), never a silent approximation.
"""

from __future__ import annotations

from typing import Any

__all__ = ["load_potts_native"]


def load_potts_native() -> Any:
    """Import and return ``ccore._potts_native``, or raise with the build hint."""
    try:
        from .ccore import _potts_native
    except ImportError as exc:  # pragma: no cover - depends on local build
        raise RuntimeError(
            "runlang='pb': the native Potts kernel (_potts_native) is not "
            "built; run `make` in src/lrgsglib/statsys/PottsModel/ccore/"
            "native/ or use runlang='py'."
        ) from exc
    return _potts_native
