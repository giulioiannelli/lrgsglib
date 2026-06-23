"""VoterModel solver backends, registered in ``statsys._solver_engine``.

Each solver is a stateless wrapper that (a) validates the model's configuration
for its backend via the model's existing ``_assert_*`` guards and (b) invokes the
corresponding kernel method. The wrappers preserve the exact call sequence the
old ``run()`` if/elif chain used, so dispatching through the registry is
behavior-preserving — the registry just becomes the single front door.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .._solver import Solver, SolverBackend
from .._solver_engine import register_solver
from .defaults import VOTER_SOLVER_NAME

if TYPE_CHECKING:
    from .VoterModel import VoterModel

__all__ = ["register_voter_solvers"]


class _VoterPySolver:
    """Pure-Python reference loop (supports the full rule + sampler grid)."""

    backend = SolverBackend.PY

    def supports(self, model: "VoterModel") -> None:
        return None

    def execute(self, model: "VoterModel", *, tqdm_on: bool = False,
                verbose: bool = False) -> None:
        model.voter_sampling(tqdm_on)

    def is_available(self) -> bool:
        return True


class _VoterPbSolver:
    """In-process pybind11 native kernel (``_voter_native``)."""

    backend = SolverBackend.PB

    def supports(self, model: "VoterModel") -> None:
        model._assert_native_supports_config("pybind")

    def execute(self, model: "VoterModel", *, tqdm_on: bool = False,
                verbose: bool = False) -> None:
        model._run_pybind()

    def is_available(self) -> bool:
        try:
            from .ccore import _voter_native  # noqa: F401
            return True
        except Exception:
            return False


class _VoterNpSolver:
    """Vectorized NumPy synchronous linear voter."""

    backend = SolverBackend.NP

    def supports(self, model: "VoterModel") -> None:
        model._assert_vectorized_supports_config("NumPy vectorized")

    def execute(self, model: "VoterModel", *, tqdm_on: bool = False,
                verbose: bool = False) -> None:
        model._run_vectorized(gpu=False)

    def is_available(self) -> bool:
        return True


class _VoterCuSolver:
    """Vectorized CuPy (GPU) synchronous linear voter."""

    backend = SolverBackend.CU

    def supports(self, model: "VoterModel") -> None:
        model._assert_vectorized_supports_config("CuPy")

    def execute(self, model: "VoterModel", *, tqdm_on: bool = False,
                verbose: bool = False) -> None:
        model._run_vectorized(gpu=True)

    def is_available(self) -> bool:
        try:
            import cupy  # noqa: F401
            return True
        except Exception:
            return False


class _VoterCSolver:
    """C subprocess backend (file-transported state)."""

    backend = SolverBackend.C

    def supports(self, model: "VoterModel") -> None:
        model._assert_native_supports_config("C subprocess")

    def execute(self, model: "VoterModel", *, tqdm_on: bool = False,
                verbose: bool = False) -> None:
        model.build_cprogram_command()
        model.run_cprogram(verbose)

    def is_available(self) -> bool:
        from ...config.lrgsg_env import LRGSG_STATSYS_VM_BIN
        return LRGSG_STATSYS_VM_BIN.exists()


_VOTER_SOLVERS: tuple[Solver, ...] = (
    _VoterPySolver(),
    _VoterPbSolver(),
    _VoterNpSolver(),
    _VoterCuSolver(),
    _VoterCSolver(),
)


def register_voter_solvers() -> None:
    """Register VoterModel's backends (idempotent; safe to call on import)."""
    for solver in _VOTER_SOLVERS:
        register_solver(VOTER_SOLVER_NAME, solver.backend, solver)


register_voter_solvers()
