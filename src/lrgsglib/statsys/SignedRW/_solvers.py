"""SignedRW solver backends, registered in ``statsys._solver_engine``.

The SignedRW package hosts two independent dynamics models, each registered
under its own name and dispatching through the shared registry (mirrors
``VoterModel/_solvers.py``):

* the single-walker family — ``py`` reference kernel + ``pb_absorb`` pybind11
  kernel (``SRW_WALKER_SOLVER_NAME``);
* the spin-copy-with-sign dynamics — ``py`` loop + ``C`` subprocess
  (``SRW_COPY_SOLVER_NAME``).

Each solver is a stateless wrapper: ``supports`` validates the configuration and
``execute`` invokes the model's existing kernel method, so dispatch is
behavior-preserving. The single native walker kernel is the absorb rule
(``pb_absorb``); that guard lives in ``SignedWalker._run_pybind`` and is
preserved verbatim.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .._solver import Solver, SolverBackend
from .._solver_engine import register_solver
from .defaults import SRW_COPY_SOLVER_NAME, SRW_WALKER_SOLVER_NAME

if TYPE_CHECKING:
    from .SignedSpinCopy import SignedSpinCopy
    from .SignedWalker import SignedWalker

__all__ = ["register_signed_rw_solvers"]


# --- single-walker family (SignedWalker / Absorbing / Killing / Sticky) ---

class _SrwWalkerPySolver:
    """Pure-Python reference walker kernel (all three rules)."""

    backend = SolverBackend.PY

    def supports(self, model: "SignedWalker") -> None:
        return None

    def execute(self, model: "SignedWalker", **kwargs: object) -> None:
        model.run_py(**kwargs)

    def is_available(self) -> bool:
        return True


class _SrwWalkerPbSolver:
    """pybind11 absorb-kernel (``runlang='pb_absorb'``, requires rule='absorb').

    The ``pb_absorb``/``rule='absorb'`` guards live in
    ``SignedWalker._run_pybind`` and fire there; ``supports`` stays a no-op so the
    native path's behaviour is preserved exactly.
    """

    backend = SolverBackend.PB

    def supports(self, model: "SignedWalker") -> None:
        return None

    def execute(self, model: "SignedWalker", **kwargs: object) -> None:
        model._run_pybind(**kwargs)

    def is_available(self) -> bool:
        from .SignedWalker import _load_srw_native
        return _load_srw_native() is not None


# --- spin-copy-with-sign dynamics (SignedSpinCopy / legacy SignedRW) -------

class _SrwCopyPySolver:
    """Pure-Python spin-copy loop."""

    backend = SolverBackend.PY

    def supports(self, model: "SignedSpinCopy") -> None:
        return None

    def execute(self, model: "SignedSpinCopy", **kwargs: object) -> None:
        model.run_py()

    def is_available(self) -> bool:
        return True


class _SrwCopyCSolver:
    """C subprocess spin-copy backend."""

    backend = SolverBackend.C

    def supports(self, model: "SignedSpinCopy") -> None:
        return None

    def execute(self, model: "SignedSpinCopy", **kwargs: object) -> None:
        model.run_c(**kwargs)

    def is_available(self) -> bool:
        return True


_WALKER_SOLVERS: tuple[Solver, ...] = (
    _SrwWalkerPySolver(),
    _SrwWalkerPbSolver(),
)
_COPY_SOLVERS: tuple[Solver, ...] = (
    _SrwCopyPySolver(),
    _SrwCopyCSolver(),
)


def register_signed_rw_solvers() -> None:
    """Register the SignedRW backends (idempotent; safe to call on import)."""
    for solver in _WALKER_SOLVERS:
        register_solver(SRW_WALKER_SOLVER_NAME, solver.backend, solver)
    for solver in _COPY_SOLVERS:
        register_solver(SRW_COPY_SOLVER_NAME, solver.backend, solver)


register_signed_rw_solvers()
