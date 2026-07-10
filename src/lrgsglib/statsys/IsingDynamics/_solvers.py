"""IsingDynamics solver backends, registered in ``statsys._solver_engine``.

Each solver is a stateless wrapper whose ``execute`` contains the *exact*
algorithm sub-dispatch (Metropolis / SA / PT / Wolff / SW / topo_met / topo_fca /
topo_cem) that the old ``run()`` if/elif chain used for the corresponding backend
family, calling the same ``_run_pybind_*`` / ``_run_cupy_*`` / sampling methods in
the same order. Dispatching through the registry is therefore behavior-preserving
-- the registry just becomes the single front door. Mirrors
``VoterModel/_solvers.py`` / ``PottsModel/_solvers.py``.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from .._solver import Solver, SolverBackend
from .._solver_engine import register_solver
from .defaults import ISING_SOLVER_NAME

if TYPE_CHECKING:
    from .IsingDynamics import IsingDynamics

__all__ = ["register_ising_solvers"]


def _resolve_modes(model: "IsingDynamics", sa_mode: bool, pt_mode: bool):
    """Recompute (rl, use_sa, use_pt) exactly as the old ``run()`` did."""
    from .IsingDynamics import _resolve_runlang

    rl = _resolve_runlang(model.runlang)
    use_sa = sa_mode or model.sa_enabled
    use_pt = pt_mode or model.pt_enabled
    return rl, use_sa, use_pt


def _is_new_style(model) -> bool:
    """True for the new scheme classes (IsingBase & leaves); the ONE registry
    key polymorphs on the instance (plan §4)."""
    from .IsingBase import IsingBase

    return isinstance(model, IsingBase)


def _reject_new_style(model, backend_name: str) -> None:
    """Hard capability error (invariant #3): the native backends for the new
    scheme classes land in Phase 2 — never a silent fallback."""
    if _is_new_style(model):
        raise NotImplementedError(
            f"The {backend_name} backend is not wired for the new Ising "
            "scheme classes yet (Phase 2: native backends); use "
            "runlang='py' or the legacy IsingDynamics."
        )


class _IsingPySolver:
    """Pure-Python reference loops (met/sa/pt/wolff/sw/topo_met/topo_fca/topo_cem)."""

    backend = SolverBackend.PY

    def supports(self, model: "IsingDynamics") -> None:
        return None

    def execute(
        self,
        model: "IsingDynamics",
        *,
        tqdm_on: bool = True,
        verbose: bool = False,
        sa_mode: bool = False,
        pt_mode: bool = False,
    ) -> None:
        if _is_new_style(model):
            # New scheme classes (IsingBase leaves): the compiled
            # ThermalEngine loop. ONE registry key, instance polymorphism.
            model._sample_py(tqdm_on=tqdm_on, verbose=verbose)
            return
        rl, use_sa, use_pt = _resolve_modes(model, sa_mode, pt_mode)
        if "topo_cem" in rl:
            model.cem_spectral_sampling(tqdm_on)
        elif "topo_met" in rl:
            model.topological_metropolis_sampling(tqdm_on)
        elif "topo_fca" in rl:
            model.topological_fca_sampling(tqdm_on)
        elif "wolff" in rl:
            model.wolff_sampling(tqdm_on)
        elif "sw" in rl:
            model.swendsen_wang_sampling(tqdm_on)
        elif use_pt:
            model.parallel_tempering_sampling(tqdm_on)
        elif use_sa:
            model.simulated_annealing_sampling(tqdm_on)
        else:
            model.metropolis_sampling(tqdm_on)

    def is_available(self) -> bool:
        return True


class _IsingPbSolver:
    """In-process pybind11 native kernel (``_ising_native``)."""

    backend = SolverBackend.PB

    def supports(self, model: "IsingDynamics") -> None:
        if _is_new_style(model):
            # New scheme classes: the leaf validates what the compiled
            # kernel can faithfully represent (hard capability errors —
            # IsingBase._pb_check_supported and overrides).
            model._pb_check_supported()
        return None

    def execute(
        self,
        model: "IsingDynamics",
        *,
        tqdm_on: bool = True,
        verbose: bool = False,
        sa_mode: bool = False,
        pt_mode: bool = False,
    ) -> None:
        if _is_new_style(model):
            # ONE registry key, instance polymorphism: the leaf owns its
            # kernel call (met/sa today; pt refused — buggy C criterion).
            model._run_pb(tqdm_on=tqdm_on, verbose=verbose)
            return
        rl, use_sa, use_pt = _resolve_modes(model, sa_mode, pt_mode)
        if "topo_cem" in rl:
            model._run_pybind_topo_cem()
        elif "topo_met" in rl:
            model._run_pybind_topo_met()
        elif "topo_fca" in rl:
            model._run_pybind_topo_fca()
        elif "wolff" in rl:
            model._run_pybind_wolff()
        elif "sw" in rl:
            model._run_pybind_sw()
        elif use_pt or rl.endswith("_pt") or "pt" in rl:
            model._run_pybind_pt()
        elif use_sa or rl.endswith("_sa") or "sa" in rl:
            model._run_pybind_sa()
        else:
            model._run_pybind_met()

    def is_available(self) -> bool:
        try:
            from .ccore import _ising_native  # noqa: F401

            return True
        except Exception:
            return False


class _IsingNpSolver:
    """Vectorized NumPy backend: sublattice-parallel sync sweeps
    (VectorSyncEngine) — new scheme classes only."""

    backend = SolverBackend.NP

    def supports(self, model: "IsingDynamics") -> None:
        if not _is_new_style(model):
            raise NotImplementedError(
                "the legacy IsingDynamics has no np backend; use the new "
                "scheme classes (IsingMetropolis, ...)."
            )
        model._np_check_supported()
        return None

    def execute(
        self,
        model: "IsingDynamics",
        *,
        tqdm_on: bool = True,
        verbose: bool = False,
        **_: object,
    ) -> None:
        model._run_np(tqdm_on=tqdm_on)

    def is_available(self) -> bool:
        return True


class _IsingCuSolver:
    """Vectorized CuPy (GPU) kernels (met/sa/pt/wolff/sw + topo fallbacks)."""

    backend = SolverBackend.CU

    def supports(self, model: "IsingDynamics") -> None:
        if _is_new_style(model):
            # Same vectorized sync engine as np, on the CuPy array module.
            model._np_check_supported()
            return None
        return None

    def execute(
        self,
        model: "IsingDynamics",
        *,
        tqdm_on: bool = True,
        verbose: bool = False,
        sa_mode: bool = False,
        pt_mode: bool = False,
    ) -> None:
        if _is_new_style(model):
            import cupy

            model._run_np(tqdm_on=tqdm_on, xp=cupy)
            return
        rl, use_sa, use_pt = _resolve_modes(model, sa_mode, pt_mode)
        if "topo_cem" in rl:
            model.cem_spectral_sampling(tqdm_on)
        elif "topo_met" in rl:
            # CuPy topological not yet implemented — fall back to Python
            model.topological_metropolis_sampling(tqdm_on)
        elif "topo_fca" in rl:
            model.topological_fca_sampling(tqdm_on)
        elif "wolff" in rl:
            model._run_cupy_wolff()
        elif "sw" in rl:
            model._run_cupy_sw()
        elif use_pt or rl.endswith("_pt") or "pt" in rl:
            model._run_cupy_pt()
        elif use_sa or rl.endswith("_sa") or "sa" in rl:
            model._run_cupy_sa()
        else:
            model._run_cupy_met()

    def is_available(self) -> bool:
        try:
            import cupy  # noqa: F401

            return True
        except Exception:
            return False


class _IsingCSolver:
    """C subprocess backend (file-transported state; met/sa/pt unified binaries)."""

    backend = SolverBackend.C

    def supports(self, model: "IsingDynamics") -> None:
        _reject_new_style(model, "C subprocess")
        return None

    def execute(
        self,
        model: "IsingDynamics",
        *,
        tqdm_on: bool = True,
        verbose: bool = False,
        sa_mode: bool = False,
        pt_mode: bool = False,
    ) -> None:
        model.build_cprogram_command()
        model.run_cprogram(verbose)

    def is_available(self) -> bool:
        bin_dir = Path(__file__).resolve().parent / "ccore" / "bin"
        return bin_dir.exists() and any(bin_dir.glob("Ising*"))


_ISING_SOLVERS: tuple[Solver, ...] = (
    _IsingPySolver(),
    _IsingPbSolver(),
    _IsingNpSolver(),
    _IsingCuSolver(),
    _IsingCSolver(),
)


def register_ising_solvers() -> None:
    """Register IsingDynamics' backends (idempotent; safe to call on import)."""
    for solver in _ISING_SOLVERS:
        register_solver(ISING_SOLVER_NAME, solver.backend, solver)


register_ising_solvers()
