"""Loader for the in-process pybind11 SIR contact-process kernel
(``_cp_native``).

One import point shared by :class:`ContactProcessSIR` and the solver
registry's availability probe. The kernel (``ccore/LRGSG_cp.c`` rate helpers
via ``ccore/cp_native.cpp``) implements ONE update -- the SIR
infection/recovery sweep of the PYTHON reference loop: every node visited
exactly once per sweep in a fresh random permutation, active nodes recover
with prob ``1 - exp(-(mu + sum_{w<0, nbr active} -w))``, inactive nodes
activate with prob ``1 - exp(-sum_{w>0, nbr active} w)``, density recorded
before each sweep. Everything the kernel cannot represent (EI dynamics,
bipolar states, per-sweep snapshots) is a hard capability error in
``_pb_check_supported`` (invariant #3), never a silent approximation.
"""

from __future__ import annotations

from typing import Any

__all__ = ["load_cp_native"]


def load_cp_native() -> Any:
    """Import and return ``ccore._cp_native``, or raise with the build hint."""
    try:
        from .ccore import _cp_native
    except ImportError as exc:  # pragma: no cover - depends on local build
        raise RuntimeError(
            "runlang='pb': the native contact-process kernel (_cp_native) is "
            "not built; run `make` in src/lrgsglib/statsys/ContactProcess/"
            "ccore/native/ or use runlang='py'."
        ) from exc
    return _cp_native
