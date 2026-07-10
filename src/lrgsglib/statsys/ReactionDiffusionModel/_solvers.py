"""ReactionDiffusionModel solver backends, registered in ``statsys._solver_engine``.

Stateless wrappers preserving the old ``run()`` call sequence so registry
dispatch is behavior-preserving. Mirrors ``KuramotoModel/_solvers.py``.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from .._solver import Solver, SolverBackend
from .._solver_engine import register_solver
from .defaults import RD_SOLVER_NAME

if TYPE_CHECKING:
    from .ReactionDiffusionModel import ReactionDiffusionModel

__all__ = ["register_rd_solvers"]


def _c_binary_exists() -> bool:
    bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    return any(bin_dir.glob("ReactionDiffusionSimulator*"))


class _RDPySolver:
    """Pure-Python RK4/Euler reference integrator."""

    backend = SolverBackend.PY

    def supports(self, model: "ReactionDiffusionModel") -> None:
        return None

    def execute(
        self,
        model: "ReactionDiffusionModel",
        *,
        tqdm_on: bool = False,
        verbose: bool = False,
    ) -> None:
        model.run_py(verbose=verbose)

    def is_available(self) -> bool:
        return True


class _RDPbSolver:
    """In-process pybind11 native RK4 kernel (``_rd_native``)."""

    backend = SolverBackend.PB

    def supports(self, model: "ReactionDiffusionModel") -> None:
        return None

    def execute(
        self,
        model: "ReactionDiffusionModel",
        *,
        tqdm_on: bool = False,
        verbose: bool = False,
    ) -> None:
        model._run_pybind()

    def is_available(self) -> bool:
        try:
            from .ccore import _rd_native  # noqa: F401

            return True
        except Exception:
            return False


class _RDCSolver:
    """C subprocess backend (file-transported state)."""

    backend = SolverBackend.C

    def supports(self, model: "ReactionDiffusionModel") -> None:
        return None

    def execute(
        self,
        model: "ReactionDiffusionModel",
        *,
        tqdm_on: bool = False,
        verbose: bool = False,
    ) -> None:
        model.build_cprogram_command()
        model.run_cprogram(verbose)

    def is_available(self) -> bool:
        return _c_binary_exists()


_RD_SOLVERS: tuple[Solver, ...] = (
    _RDPySolver(),
    _RDPbSolver(),
    _RDCSolver(),
)


def register_rd_solvers() -> None:
    """Register ReactionDiffusionModel's backends (idempotent; safe on import)."""
    for solver in _RD_SOLVERS:
        register_solver(RD_SOLVER_NAME, solver.backend, solver)


register_rd_solvers()
