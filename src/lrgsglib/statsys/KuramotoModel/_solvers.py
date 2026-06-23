"""KuramotoModel solver backends, registered in ``statsys._solver_engine``.

Stateless wrappers preserving the old ``run()`` call sequence so registry
dispatch is behavior-preserving. Mirrors ``VoterModel/_solvers.py``.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from .._solver import Solver, SolverBackend
from .._solver_engine import register_solver
from .defaults import KURAMOTO_SOLVER_NAME

if TYPE_CHECKING:
    from .KuramotoModel import KuramotoModel

__all__ = ["register_kuramoto_solvers"]


def _c_binary_exists() -> bool:
    bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    return any(bin_dir.glob("KuramotoSimulator*"))


class _KuramotoPySolver:
    """Pure-Python RK4/Euler reference integrator."""

    backend = SolverBackend.PY

    def supports(self, model: "KuramotoModel") -> None:
        return None

    def execute(self, model: "KuramotoModel", *, verbose: bool = False) -> None:
        model.run_py(verbose=verbose)

    def is_available(self) -> bool:
        return True


class _KuramotoPbSolver:
    """In-process pybind11 native RK4 kernel (``_kuramoto_native``)."""

    backend = SolverBackend.PB

    def supports(self, model: "KuramotoModel") -> None:
        return None

    def execute(self, model: "KuramotoModel", *, verbose: bool = False) -> None:
        model._run_pybind()

    def is_available(self) -> bool:
        try:
            from .ccore import _kuramoto_native  # noqa: F401
            return True
        except Exception:
            return False


class _KuramotoCSolver:
    """C subprocess backend (file-transported state)."""

    backend = SolverBackend.C

    def supports(self, model: "KuramotoModel") -> None:
        return None

    def execute(self, model: "KuramotoModel", *, verbose: bool = False) -> None:
        model.build_cprogram_command()
        model.run_cprogram(verbose)

    def is_available(self) -> bool:
        return _c_binary_exists()


_KURAMOTO_SOLVERS: tuple[Solver, ...] = (
    _KuramotoPySolver(),
    _KuramotoPbSolver(),
    _KuramotoCSolver(),
)


def register_kuramoto_solvers() -> None:
    """Register KuramotoModel's backends (idempotent; safe to call on import)."""
    for solver in _KURAMOTO_SOLVERS:
        register_solver(KURAMOTO_SOLVER_NAME, solver.backend, solver)


register_kuramoto_solvers()
