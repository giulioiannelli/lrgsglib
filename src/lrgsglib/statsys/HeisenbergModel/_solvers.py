"""HeisenbergModel solver backends, registered in ``statsys._solver_engine``.

Each solver is a stateless wrapper that invokes the corresponding kernel method,
preserving the exact call sequence the old ``run()`` if/elif chain used so that
dispatching through the registry is behavior-preserving -- the registry just
becomes the single front door. Mirrors ``XYModel/_solvers.py``.

ONE registry key ("heisenberg") serves BOTH the legacy ``HeisenbergModel`` and
the new scheme classes (``HeisenbergBase`` leaves): the solver polymorphs on
the instance (statsys unification D-B4), exactly as ``XYModel/_solvers.py``
does. The C-subprocess backend stays legacy-only (decision D-B5, frozen argv).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from .._solver import Solver, SolverBackend
from .._solver_engine import register_solver
from .defaults import HEISENBERG_SOLVER_NAME

if TYPE_CHECKING:
    from .HeisenbergModel import HeisenbergModel

__all__ = ["register_heisenberg_solvers"]


def _is_new_style(model) -> bool:
    """True for the new scheme classes (HeisenbergBase & leaves); the ONE
    registry key polymorphs on the instance (D-B4)."""
    from .HeisenbergBase import HeisenbergBase

    return isinstance(model, HeisenbergBase)


class _HeisenbergPySolver:
    """Pure-Python Metropolis reference loop."""

    backend = SolverBackend.PY

    def supports(self, model: "HeisenbergModel") -> None:
        return None

    def execute(
        self,
        model: "HeisenbergModel",
        *,
        tqdm_on: bool = False,
        verbose: bool = False,
    ) -> None:
        if _is_new_style(model):
            # New scheme classes (HeisenbergBase leaves): the compiled shared
            # ThermalEngine loop. ONE registry key, instance polymorphism.
            model._sample_py(tqdm_on=tqdm_on, verbose=verbose)
            return
        model.run_py(verbose=verbose)

    def is_available(self) -> bool:
        return True


class _HeisenbergPbSolver:
    """In-process pybind11 native kernel (``_heisenberg_native``)."""

    backend = SolverBackend.PB

    def supports(self, model: "HeisenbergModel") -> None:
        if _is_new_style(model):
            # New scheme classes: the leaf validates what the compiled
            # kernel can faithfully represent (hard capability errors).
            model._pb_check_supported()
        return None

    def execute(
        self,
        model: "HeisenbergModel",
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
            from .ccore import _heisenberg_native  # noqa: F401

            return True
        except Exception:
            return False


class _HeisenbergNpSolver:
    """Vectorized NumPy backend: sublattice-parallel sync sweeps
    (VectorSyncEngine) — new scheme classes only."""

    backend = SolverBackend.NP

    def supports(self, model: "HeisenbergModel") -> None:
        if not _is_new_style(model):
            raise NotImplementedError(
                "the legacy HeisenbergModel has no np backend; use the new "
                "scheme classes (HeisenbergMetropolis, ...)."
            )
        model._np_check_supported()
        return None

    def execute(
        self,
        model: "HeisenbergModel",
        *,
        tqdm_on: bool = False,
        verbose: bool = False,
        **_: object,
    ) -> None:
        model._run_np(tqdm_on=tqdm_on)

    def is_available(self) -> bool:
        return True


class _HeisenbergCuSolver:
    """Vectorized CuPy (GPU) backend: the same VectorSyncEngine as np,
    on the CuPy array module — new scheme classes only."""

    backend = SolverBackend.CU

    def supports(self, model: "HeisenbergModel") -> None:
        if not _is_new_style(model):
            raise NotImplementedError(
                "the legacy HeisenbergModel has no cu backend; use the new "
                "scheme classes (HeisenbergMetropolis, ...)."
            )
        model._np_check_supported()
        return None

    def execute(
        self,
        model: "HeisenbergModel",
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


class _HeisenbergCSolver:
    """C subprocess backend (file-transported state)."""

    backend = SolverBackend.C

    def supports(self, model: "HeisenbergModel") -> None:
        if _is_new_style(model):
            # Hard capability error (invariant #3): the C subprocess stays
            # on the untouched legacy class (D-B5) — never a silent fallback.
            raise NotImplementedError(
                "The C-subprocess backend is not wired for the new "
                "Heisenberg scheme classes; use runlang='py'/'pb' or the "
                "legacy HeisenbergModel."
            )
        return None

    def execute(
        self,
        model: "HeisenbergModel",
        *,
        tqdm_on: bool = False,
        verbose: bool = False,
    ) -> None:
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
    _HeisenbergNpSolver(),
    _HeisenbergCuSolver(),
    _HeisenbergCSolver(),
)


def register_heisenberg_solvers() -> None:
    """Register HeisenbergModel's backends (idempotent; safe on import)."""
    for solver in _HEISENBERG_SOLVERS:
        register_solver(HEISENBERG_SOLVER_NAME, solver.backend, solver)


register_heisenberg_solvers()
