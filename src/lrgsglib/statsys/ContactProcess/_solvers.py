"""ContactProcess solver backends, registered in ``statsys._solver_engine``.

Mirrors ``VoterModel/_solvers.py``: each solver is a stateless wrapper that
validates the model's configuration (``supports``) and invokes its existing
kernel path (``execute``), so dispatching through the registry is
behavior-preserving — the registry just becomes the single front door.

ContactProcess exposes three families: the pure-Python loop (selected
polymorphically via ``run_py`` — the scalar ``ds1step`` sweep for SIR, the
numba lambda-cache sweep for EI), the in-process pybind11 kernel
(``pb``, SIR only — EI is a hard capability error), and the C subprocess
(``C0`` SIR / ``C1*`` EI).
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from .._solver import Solver, SolverBackend
from .._solver_engine import register_solver
from .defaults import CP_SOLVER_NAME

if TYPE_CHECKING:
    from .ContactProcess import ContactProcessBase

__all__ = ["register_contact_process_solvers"]


class _CpPySolver:
    """Pure-Python loop (SIR scalar sweep / EI numba lambda-cache, selected
    polymorphically via ``_sample_py``)."""

    backend = SolverBackend.PY

    def supports(self, model: "ContactProcessBase") -> None:
        return None

    def execute(
        self,
        model: "ContactProcessBase",
        *,
        tqdm_on: bool = False,
        verbose: bool = False,
    ) -> None:
        # Bare sampling loop: the run() front door (RunHostMixin) owns the
        # output lifecycle around the solver.
        model._sample_py(tqdm_on=tqdm_on)

    def is_available(self) -> bool:
        return True


class _CpPbSolver:
    """In-process pybind11 native kernel (``_cp_native``; SIR only — the
    model's ``_pb_check_supported`` raises the hard NotImplementedError
    for everything the kernel cannot represent, EI included)."""

    backend = SolverBackend.PB

    def supports(self, model: "ContactProcessBase") -> None:
        model._pb_check_supported()

    def execute(
        self,
        model: "ContactProcessBase",
        *,
        tqdm_on: bool = False,
        verbose: bool = False,
    ) -> None:
        model._run_pb(tqdm_on=tqdm_on, verbose=verbose)

    def is_available(self) -> bool:
        try:
            from .ccore import _cp_native  # noqa: F401

            return True
        except Exception:
            return False


class _CpCSolver:
    """C subprocess backend (file-transported state; SIR ``C0`` / EI ``C1*``)."""

    backend = SolverBackend.C

    def supports(self, model: "ContactProcessBase") -> None:
        return None

    def execute(
        self,
        model: "ContactProcessBase",
        *,
        tqdm_on: bool = False,
        verbose: bool = False,
    ) -> None:
        model.build_cprogram_command()
        model.run_cprogram(verbose)

    def is_available(self) -> bool:
        bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
        return bin_dir.exists() and any(bin_dir.glob("ContactProcess*"))


_CP_SOLVERS: tuple[Solver, ...] = (
    _CpPySolver(),
    _CpPbSolver(),
    _CpCSolver(),
)


def register_contact_process_solvers() -> None:
    """Register ContactProcess backends (idempotent; safe to call on import)."""
    for solver in _CP_SOLVERS:
        register_solver(CP_SOLVER_NAME, solver.backend, solver)


register_contact_process_solvers()
