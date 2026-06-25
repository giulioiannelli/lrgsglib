"""XYModel solver backends, registered in ``statsys._solver_engine``.

Each solver is a stateless wrapper that invokes the corresponding kernel method,
preserving the exact call sequence the old ``run()`` if/elif chain used so that
dispatching through the registry is behavior-preserving -- the registry just
becomes the single front door. Mirrors ``PottsModel/_solvers.py``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .._solver import Solver, SolverBackend
from .._solver_engine import register_solver
from .defaults import XY_SOLVER_NAME

if TYPE_CHECKING:
    from .XYModel import XYModel

__all__ = ["register_xy_solvers"]


class _XYPySolver:
    """Pure-Python Metropolis reference loop."""

    backend = SolverBackend.PY

    def supports(self, model: "XYModel") -> None:
        return None

    def execute(self, model: "XYModel", *, verbose: bool = False) -> None:
        model.run_py(verbose=verbose)

    def is_available(self) -> bool:
        return True


class _XYPbSolver:
    """In-process pybind11 native kernel (``_xy_native``)."""

    backend = SolverBackend.PB

    def supports(self, model: "XYModel") -> None:
        return None

    def execute(self, model: "XYModel", *, verbose: bool = False) -> None:
        model._run_pybind()

    def is_available(self) -> bool:
        try:
            from .ccore import _xy_native  # noqa: F401
            return True
        except Exception:
            return False


class _XYCSolver:
    """C subprocess backend (file-transported state)."""

    backend = SolverBackend.C

    def supports(self, model: "XYModel") -> None:
        return None

    def execute(self, model: "XYModel", *, verbose: bool = False) -> None:
        model.build_cprogram_command()
        model.run_cprogram(verbose)

    def is_available(self) -> bool:
        return model_c_binary_exists()


def model_c_binary_exists() -> bool:
    """Whether the compiled XYSimulator C binary is present."""
    from pathlib import Path

    bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    return any(bin_dir.glob("XYSimulator*"))


_XY_SOLVERS: tuple[Solver, ...] = (
    _XYPySolver(),
    _XYPbSolver(),
    _XYCSolver(),
)


def register_xy_solvers() -> None:
    """Register XYModel's backends (py/pb/C) in the shared solver registry."""
    for solver in _XY_SOLVERS:
        register_solver(XY_SOLVER_NAME, solver.backend, solver)


register_xy_solvers()
