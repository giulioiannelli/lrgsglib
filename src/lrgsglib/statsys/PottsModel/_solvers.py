"""PottsModel solver backends, registered in ``statsys._solver_engine``.

Each solver is a stateless wrapper that invokes the corresponding kernel method,
preserving the exact call sequence the old ``run()`` if/elif chain used so that
dispatching through the registry is behavior-preserving -- the registry just
becomes the single front door. Mirrors ``VoterModel/_solvers.py``.

ONE registry key ("potts") serves BOTH the legacy ``PottsModel`` and the new
scheme classes (``PottsBase`` leaves): the solver polymorphs on the instance
(statsys unification D-B4), exactly as ``IsingDynamics/_solvers.py`` does.
The C-subprocess backend stays legacy-only (decision D-B5, frozen argv).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .._solver import Solver, SolverBackend
from .._solver_engine import register_solver
from .defaults import POTTS_SOLVER_NAME

if TYPE_CHECKING:
    from .PottsModel import PottsModel

__all__ = ["register_potts_solvers"]


def _is_new_style(model) -> bool:
    """True for the new scheme classes (PottsBase & leaves); the ONE registry
    key polymorphs on the instance (D-B4)."""
    from .PottsBase import PottsBase

    return isinstance(model, PottsBase)


class _PottsPySolver:
    """Pure-Python Metropolis reference loop."""

    backend = SolverBackend.PY

    def supports(self, model: "PottsModel") -> None:
        return None

    def execute(
        self,
        model: "PottsModel",
        *,
        tqdm_on: bool = False,
        verbose: bool = False,
    ) -> None:
        if _is_new_style(model):
            # New scheme classes (PottsBase leaves): the compiled shared
            # ThermalEngine loop. ONE registry key, instance polymorphism.
            model._sample_py(tqdm_on=tqdm_on, verbose=verbose)
            return
        model.run_py(verbose=verbose)

    def is_available(self) -> bool:
        return True


class _PottsPbSolver:
    """In-process pybind11 native kernel (``_potts_native``)."""

    backend = SolverBackend.PB

    def supports(self, model: "PottsModel") -> None:
        if _is_new_style(model):
            # New scheme classes: the leaf validates what the compiled
            # kernel can faithfully represent (hard capability errors).
            model._pb_check_supported()
        return None

    def execute(
        self,
        model: "PottsModel",
        *,
        tqdm_on: bool = False,
        verbose: bool = False,
    ) -> None:
        if _is_new_style(model):
            model._run_pb(tqdm_on=tqdm_on, verbose=verbose)
            return
        model._run_pybind()

    def is_available(self) -> bool:
        try:
            from .ccore import _potts_native  # noqa: F401

            return True
        except Exception:
            return False


class _PottsNpSolver:
    """Vectorized NumPy backend: sublattice-parallel sync sweeps
    (VectorSyncEngine) — new scheme classes only."""

    backend = SolverBackend.NP

    def supports(self, model: "PottsModel") -> None:
        if not _is_new_style(model):
            raise NotImplementedError(
                "the legacy PottsModel has no np backend; use the new "
                "scheme classes (PottsMetropolis, ...)."
            )
        model._np_check_supported()
        return None

    def execute(
        self,
        model: "PottsModel",
        *,
        tqdm_on: bool = False,
        verbose: bool = False,
        **_: object,
    ) -> None:
        model._run_np(tqdm_on=tqdm_on)

    def is_available(self) -> bool:
        return True


class _PottsCuSolver:
    """Vectorized CuPy (GPU) backend: the same VectorSyncEngine as np,
    on the CuPy array module — new scheme classes only."""

    backend = SolverBackend.CU

    def supports(self, model: "PottsModel") -> None:
        if not _is_new_style(model):
            raise NotImplementedError(
                "the legacy PottsModel has no cu backend; use the new "
                "scheme classes (PottsMetropolis, ...)."
            )
        model._np_check_supported()
        return None

    def execute(
        self,
        model: "PottsModel",
        *,
        tqdm_on: bool = False,
        verbose: bool = False,
        **_: object,
    ) -> None:
        import cupy

        model._run_np(tqdm_on=tqdm_on, xp=cupy)

    def is_available(self) -> bool:
        try:
            import cupy  # noqa: F401

            return True
        except Exception:
            return False


class _PottsCSolver:
    """C subprocess backend (file-transported state)."""

    backend = SolverBackend.C

    def supports(self, model: "PottsModel") -> None:
        if _is_new_style(model):
            # Hard capability error (invariant #3): the C subprocess stays
            # on the untouched legacy class (D-B5) — never a silent fallback.
            raise NotImplementedError(
                "The C-subprocess backend is not wired for the new Potts "
                "scheme classes; use runlang='py'/'pb' or the legacy "
                "PottsModel."
            )
        return None

    def execute(
        self,
        model: "PottsModel",
        *,
        tqdm_on: bool = False,
        verbose: bool = False,
    ) -> None:
        model.build_cprogram_command()
        model.run_cprogram(verbose)

    def is_available(self) -> bool:
        return model_c_binary_exists()


def model_c_binary_exists() -> bool:
    """Whether the compiled PottsSimulator C binary is present."""
    from pathlib import Path

    bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    return any(bin_dir.glob("PottsSimulator*"))


_POTTS_SOLVERS: tuple[Solver, ...] = (
    _PottsPySolver(),
    _PottsPbSolver(),
    _PottsNpSolver(),
    _PottsCuSolver(),
    _PottsCSolver(),
)


def register_potts_solvers() -> None:
    """Register PottsModel's backends (idempotent; safe to call on import)."""
    for solver in _POTTS_SOLVERS:
        register_solver(POTTS_SOLVER_NAME, solver.backend, solver)


register_potts_solvers()
