"""HeisenbergModel solver backends, registered in ``statsys._solver_engine``.

Each solver is a stateless wrapper that invokes the corresponding kernel method,
preserving the exact call sequence the old ``run()`` if/elif chain used so that
dispatching through the registry is behavior-preserving -- the registry just
becomes the single front door. Mirrors ``PottsModel/_solvers.py``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .._solver import Solver, SolverBackend
from .._solver_engine import register_solver
from .defaults import HEISENBERG_SOLVER_NAME

if TYPE_CHECKING:
    from .HeisenbergModel import HeisenbergModel

__all__ = ["register_heisenberg_solvers"]


class _HeisenbergPySolver:
    """Pure-Python Metropolis reference loop."""

    backend = SolverBackend.PY

    def supports(self, model: "HeisenbergModel") -> None:
        return None

    def execute(self, model: "HeisenbergModel", *, verbose: bool = False) -> None:
        model.run_py(verbose=verbose)

    def is_available(self) -> bool:
        return True


class _HeisenbergPbSolver:
    """In-process pybind11 native kernel (``_heisenberg_native``)."""

    backend = SolverBackend.PB

    def supports(self, model: "HeisenbergModel") -> None:
        return None

    def execute(self, model: "HeisenbergModel", *, verbose: bool = False) -> None:
        model._run_pybind()

    def is_available(self) -> bool:
        try:
            from .ccore import _heisenberg_native  # noqa: F401
            return True
        except Exception:
            return False


class _HeisenbergCSolver:
    """C subprocess backend (file-transported state)."""

    backend = SolverBackend.C

    def supports(self, model: "HeisenbergModel") -> None:
        return None

    def execute(self, model: "HeisenbergModel", *, verbose: bool = False) -> None:
        model.build_cprogram_command()
        model.run_cprogram(verbose)

    def is_available(self) -> bool:
        return model_c_binary_exists()


def model_c_binary_exists() -> bool:
    """Whether the compiled HeisenbergSimulator C binary is present."""
    from pathlib import Path

    bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    return any(bin_dir.glob("HeisenbergSimulator*"))


_HEISENBERG_SOLVERS: tuple[Solver, ...] = (
    _HeisenbergPySolver(),
    _HeisenbergPbSolver(),
    _HeisenbergCSolver(),
)


def register_heisenberg_solvers() -> None:
    """Register HeisenbergModel's backends (idempotent; safe to call on import)."""
    for solver in _HEISENBERG_SOLVERS:
        register_solver(HEISENBERG_SOLVER_NAME, solver.backend, solver)


register_heisenberg_solvers()
