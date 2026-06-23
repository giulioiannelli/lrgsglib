"""CoupledODEModel solver backends, registered in ``statsys._solver_engine``.

Stateless wrappers preserving the old ``run()`` call sequence so registry
dispatch is behavior-preserving. Mirrors ``KuramotoModel/_solvers.py``.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from .._solver import Solver, SolverBackend
from .._solver_engine import register_solver
from .defaults import CODE_SOLVER_NAME

if TYPE_CHECKING:
    from .CoupledODEModel import CoupledODEModel

__all__ = ["register_code_solvers"]


def _c_binary_exists() -> bool:
    bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    return any(bin_dir.glob("CoupledODESimulator*"))


class _CodePySolver:
    """Pure-Python RK4/Euler reference integrator."""

    backend = SolverBackend.PY

    def supports(self, model: "CoupledODEModel") -> None:
        return None

    def execute(self, model: "CoupledODEModel", *, verbose: bool = False) -> None:
        model.run_py(verbose=verbose)

    def is_available(self) -> bool:
        return True


class _CodePbSolver:
    """In-process pybind11 native RK4 kernel (``_code_native``)."""

    backend = SolverBackend.PB

    def supports(self, model: "CoupledODEModel") -> None:
        return None

    def execute(self, model: "CoupledODEModel", *, verbose: bool = False) -> None:
        model._run_pybind()

    def is_available(self) -> bool:
        try:
            from .ccore import _code_native  # noqa: F401
            return True
        except Exception:
            return False


class _CodeCSolver:
    """C subprocess backend (file-transported state)."""

    backend = SolverBackend.C

    def supports(self, model: "CoupledODEModel") -> None:
        return None

    def execute(self, model: "CoupledODEModel", *, verbose: bool = False) -> None:
        model.build_cprogram_command()
        model.run_cprogram(verbose)

    def is_available(self) -> bool:
        return _c_binary_exists()


_CODE_SOLVERS: tuple[Solver, ...] = (
    _CodePySolver(),
    _CodePbSolver(),
    _CodeCSolver(),
)


def register_code_solvers() -> None:
    """Register CoupledODEModel's backends (idempotent; safe to call on import)."""
    for solver in _CODE_SOLVERS:
        register_solver(CODE_SOLVER_NAME, solver.backend, solver)


register_code_solvers()
