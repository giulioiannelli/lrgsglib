"""MultiSpeciesModel solver backends, registered in ``statsys._solver_engine``.

Each solver is a stateless wrapper that invokes the corresponding kernel method,
preserving the exact call sequence the old ``run()`` if/elif chain used so that
dispatching through the registry is behavior-preserving -- the registry just
becomes the single front door. Mirrors ``PottsModel/_solvers.py``.

ONE registry key ("multispec") serves BOTH the legacy ``MultiSpeciesModel``
and the new scheme classes (``MultiSpeciesBase`` leaves): the solver
polymorphs on the instance (statsys unification D-B4), exactly as
``PottsModel/_solvers.py`` does. The C-subprocess backend stays legacy-only
(decision D-B5, frozen argv).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .._solver import Solver, SolverBackend
from .._solver_engine import register_solver
from .defaults import MULTISPEC_SOLVER_NAME

if TYPE_CHECKING:
    from .MultiSpeciesModel import MultiSpeciesModel

__all__ = ["register_multispec_solvers"]


def _is_new_style(model) -> bool:
    """True for the new scheme classes (MultiSpeciesBase & leaves); the ONE
    registry key polymorphs on the instance (D-B4)."""
    from .MultiSpeciesBase import MultiSpeciesBase

    return isinstance(model, MultiSpeciesBase)


class _MultiSpecPySolver:
    """Pure-Python Metropolis reference loop."""

    backend = SolverBackend.PY

    def supports(self, model: "MultiSpeciesModel") -> None:
        return None

    def execute(
        self,
        model: "MultiSpeciesModel",
        *,
        tqdm_on: bool = False,
        verbose: bool = False,
    ) -> None:
        if _is_new_style(model):
            # New scheme classes (MultiSpeciesBase leaves): the compiled
            # shared ThermalEngine loop. ONE registry key, instance
            # polymorphism.
            model._sample_py(tqdm_on=tqdm_on, verbose=verbose)
            return
        model.run_py(verbose=verbose)

    def is_available(self) -> bool:
        return True


class _MultiSpecPbSolver:
    """In-process pybind11 native kernel (``_multispec_native``)."""

    backend = SolverBackend.PB

    def supports(self, model: "MultiSpeciesModel") -> None:
        if _is_new_style(model):
            # New scheme classes: the leaf validates what the compiled
            # kernel can faithfully represent (hard capability errors).
            model._pb_check_supported()
        return None

    def execute(
        self,
        model: "MultiSpeciesModel",
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
            from .ccore import _multispec_native  # noqa: F401

            return True
        except Exception:
            return False


class _MultiSpecCSolver:
    """C subprocess backend (file-transported state)."""

    backend = SolverBackend.C

    def supports(self, model: "MultiSpeciesModel") -> None:
        if _is_new_style(model):
            # Hard capability error (invariant #3): the C subprocess stays
            # on the untouched legacy class (D-B5) — never a silent fallback.
            raise NotImplementedError(
                "The C-subprocess backend is not wired for the new "
                "multi-species scheme classes; use runlang='py'/'pb' or the "
                "legacy MultiSpeciesModel."
            )
        return None

    def execute(
        self,
        model: "MultiSpeciesModel",
        *,
        tqdm_on: bool = False,
        verbose: bool = False,
    ) -> None:
        model.build_cprogram_command()
        model.run_cprogram(verbose)

    def is_available(self) -> bool:
        return model_c_binary_exists()


def model_c_binary_exists() -> bool:
    """Whether the compiled MultiSpeciesSimulator C binary is present."""
    from pathlib import Path

    bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
    return any(bin_dir.glob("MultiSpeciesSimulator*"))


_MULTISPEC_SOLVERS: tuple[Solver, ...] = (
    _MultiSpecPySolver(),
    _MultiSpecPbSolver(),
    _MultiSpecCSolver(),
)


def register_multispec_solvers() -> None:
    """Register MultiSpeciesModel's backends (idempotent; safe on import)."""
    for solver in _MULTISPEC_SOLVERS:
        register_solver(MULTISPEC_SOLVER_NAME, solver.backend, solver)


register_multispec_solvers()
