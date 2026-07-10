"""IsingSimulatedAnnealing — cooling-schedule optimization (Layer-1 scheme).

One system driven through a decreasing temperature schedule by the shared
``statsys._thermal.ThermalEngine``: the engine is recompiled once per
temperature stage (configure-then-run, plan §5 — the sweep stays branch-free
inside a stage). Acceptance ``rule`` and async ``order`` are the same Layer-2
axes as on :class:`IsingMetropolis` (citations in
``.agents/ref/00-REFERENCES.md`` §MCMC).

The legacy TFCA algorithm ("topological field-cooling annealing") is this
scheme with ``field='spectral'`` (plan §10): the external field is built at
init from the softmax-weighted low-RBIM-energy eigenvectors
(``statsys._spectral.build_spectral_field``) and enters ΔE and E first-class,
so the quench relaxes toward the spectral sign structure.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np
import tqdm

from .._spectral import (
    build_spectral_field,
    eigvec_subspace_matrix,
    rbim_subspace_energies,
)
from .._thermal import THERMAL_RULES, ThermalEngine, resolve_tie_flip
from ._native import PB_ORDER_MAP, load_ising_native
from .defaults import (
    ISING_FIELD_MODE_DEFAULT,
    ISING_FIELD_MODES,
    ISING_ORDER_DEFAULT,
    ISING_RULE_DEFAULT,
    ISING_SA_COOLING_RATE_DEFAULT,
    ISING_SA_N_TEMPERATURES_DEFAULT,
    ISING_SA_SCHEDULE_DEFAULT,
    ISING_SA_SCHEDULES,
    ISING_SA_STEPS_PER_T_DEFAULT,
    ISING_SA_T_FINAL_DEFAULT,
    ISING_SA_T_INIT_DEFAULT,
    ISING_SPECTRAL_FIELD_STRENGTH_DEFAULT,
    ISING_SPECTRAL_N_MODES_DEFAULT,
    ISING_SPECTRAL_TAU_DEFAULT,
    ISING_TIE_FLIP_DEFAULT,
    ISING_UPD_MODE_DEFAULT,
)
from .IsingBase import IsingBase

if TYPE_CHECKING:
    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph

__all__ = ["IsingSimulatedAnnealing", "generate_temperature_schedule"]


def generate_temperature_schedule(
    schedule: str,
    T_init: float,
    T_final: float,
    n_temperatures: int,
    cooling_rate: float,
) -> np.ndarray:
    """Build the stage-temperature array for a named cooling schedule.

    'linear' interpolates ``T_init → T_final``; 'exponential' is
    ``T_init · rate^k`` (``T_final`` unused — legacy semantics preserved);
    'logarithmic' is ``T_init / ln(2 + k)``.
    """
    if schedule not in ISING_SA_SCHEDULES:
        raise ValueError(
            f"schedule must be one of {ISING_SA_SCHEDULES}; got {schedule!r}."
        )
    n = int(n_temperatures)
    if n < 1:
        raise ValueError(f"n_temperatures must be >= 1; got {n}.")
    if schedule == "linear":
        return np.linspace(T_init, T_final, n)
    if schedule == "exponential":
        return T_init * (cooling_rate ** np.arange(n))
    return T_init / np.log(2 + np.arange(n))  # logarithmic


class IsingSimulatedAnnealing(IsingBase):
    """Simulated annealing on the shared MCMC engine.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph carrying the couplings.
    T_init, T_final : float
        Endpoints of the cooling schedule (``k_B = 1``).
    schedule : {'linear', 'exponential', 'logarithmic'}
        Named cooling schedule (see :func:`generate_temperature_schedule`).
    cooling_rate : float
        Per-stage factor of the exponential schedule.
    n_temperatures : int
        Number of temperature stages.
    steps_per_T : int
        Monte Carlo sweeps per stage; the total step budget is
        ``n_temperatures * steps_per_T`` (fixed by the schedule — ``run()``
        rejects a ``steps=`` override).
    T_schedule : ndarray, optional
        Explicit stage temperatures; overrides the named schedule.
    rule, order, tie_flip_p
        Same Layer-2 kinetics axes as :class:`IsingMetropolis`.
    field : {'uniform', 'spectral'} or ndarray
        External-field construction: an array is used as-is (the 'uniform'
        mode), ``'spectral'`` builds the softmax-weighted eigenvector field
        at init (the legacy TFCA field). Default: zero field.
    spectral_n_modes, spectral_tau, spectral_field_strength
        Parameters of the spectral field (``field='spectral'`` only).
    **kwargs
        Forwarded to :class:`IsingBase` (``observables``, ``coupling_norm``,
        ``ic``, ``seed``, ``savedisk``, ...).

    After a run, the per-stage curve is available as :attr:`sa_temps`,
    :attr:`sa_energy` (intensive) and :attr:`sa_magn` (see :meth:`sa_curve`);
    the sweep-resolution traces are the standard observables.
    """

    def __init__(
        self,
        sg: "SignedGraph",
        *,
        T_init: float = ISING_SA_T_INIT_DEFAULT,
        T_final: float = ISING_SA_T_FINAL_DEFAULT,
        schedule: str = ISING_SA_SCHEDULE_DEFAULT,
        cooling_rate: float = ISING_SA_COOLING_RATE_DEFAULT,
        n_temperatures: int = ISING_SA_N_TEMPERATURES_DEFAULT,
        steps_per_T: int = ISING_SA_STEPS_PER_T_DEFAULT,
        T_schedule: np.ndarray | None = None,
        rule: str = ISING_RULE_DEFAULT,
        order: str = ISING_ORDER_DEFAULT,
        tie_flip_p: float | str | None = None,
        field: str | np.ndarray | None = ISING_FIELD_MODE_DEFAULT,
        spectral_n_modes: int = ISING_SPECTRAL_N_MODES_DEFAULT,
        spectral_tau: float = ISING_SPECTRAL_TAU_DEFAULT,
        spectral_field_strength: float = ISING_SPECTRAL_FIELD_STRENGTH_DEFAULT,
        **kwargs: Any,
    ) -> None:
        if rule not in THERMAL_RULES:
            raise ValueError(
                f"rule must be one of {THERMAL_RULES}; got {rule!r}."
            )
        if T_schedule is not None:
            stage_T = np.asarray(T_schedule, dtype=np.float64)
        else:
            stage_T = generate_temperature_schedule(
                schedule, T_init, T_final, n_temperatures, cooling_rate
            )
        if stage_T.ndim != 1 or len(stage_T) < 1:
            raise ValueError("T_schedule must be a non-empty 1-D array.")
        if np.any(stage_T < 0.0):
            raise ValueError(
                "All schedule temperatures must be >= 0; got "
                f"min={stage_T.min()}."
            )

        # The field axis: a string selects the construction mode, an array
        # IS the uniform mode (passed through to the DynSys field slot).
        if isinstance(field, str):
            if field not in ISING_FIELD_MODES:
                raise ValueError(
                    f"field must be one of {ISING_FIELD_MODES} or an "
                    f"explicit array; got {field!r}."
                )
            field_mode, field_arr = field, None
        else:
            field_mode, field_arr = "uniform", field

        self.T_schedule = stage_T
        # The schedule's identity (name + generator parameters), kept for
        # the native backend: the C kernel regenerates the temperature grid
        # itself, so the pb guard must know which named schedule this is.
        self.schedule = "custom" if T_schedule is not None else schedule
        self.T_init = float(T_init)
        self.T_final = float(T_final)
        self.cooling_rate = float(cooling_rate)
        self.n_temperatures = len(stage_T)
        self.steps_per_T = int(steps_per_T)
        if self.steps_per_T < 1:
            raise ValueError(f"steps_per_T must be >= 1; got {steps_per_T}.")
        self.T = float(stage_T[0])
        self.rule = rule
        self.upd_mode = ISING_UPD_MODE_DEFAULT
        self.order = order
        self.tie_flip_p = resolve_tie_flip(
            rule, tie_flip_p, default=ISING_TIE_FLIP_DEFAULT
        )
        self.field_mode = field_mode
        self.spectral_n_modes = int(spectral_n_modes)
        self.spectral_tau = float(spectral_tau)
        self.spectral_field_strength = float(spectral_field_strength)
        self.spectral_field: np.ndarray | None = None
        self.sa_temps: np.ndarray | None = None
        self.sa_energy: np.ndarray | None = None
        self.sa_magn: np.ndarray | None = None

        super().__init__(
            sg,
            field=field_arr,
            steps=self.n_temperatures * self.steps_per_T,
            **kwargs,
        )

        # Setup-time capability check (engine construction validates the
        # rule/order/tie combination NOW, not mid-run).
        self._make_engine()

    # ------------------------------------------------------------------
    # Run-dirname schema (Phase C)
    # ------------------------------------------------------------------
    _sch_token = "sa"

    def _ns_token_value(self) -> int | None:
        return None  # horizon derived from the scheme numerics

    def _physics_tokens(self) -> list:
        return [
            ("Ti", float(self.T_init)),
            ("Tf", float(self.T_final)),
            ("nT", int(self.n_temperatures)),
            ("spt", int(self.steps_per_T)),
            self._field_token(),
        ]

    def _axis_tokens(self) -> list:
        return [("sched", self.schedule, "exponential")]

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------
    def init_dynamics(self, custom: Any = None) -> None:
        super().init_dynamics(custom)
        if self.field_mode == "spectral":
            n_modes = min(self.spectral_n_modes, self.N)
            V = eigvec_subspace_matrix(self.sg, n_modes)
            rbim_E = rbim_subspace_energies(self.sg, n_modes)
            self.spectral_field = build_spectral_field(
                V, rbim_E, self.spectral_tau, self.spectral_field_strength
            )
            self.field = self.spectral_field

    def initialize_run_parameters(
        self, steps: int | None = None, simref: float | None = None
    ) -> None:
        if steps is not None or simref is not None:
            raise ValueError(
                "IsingSimulatedAnnealing derives its step budget from the "
                "schedule (n_temperatures x steps_per_T); configure those "
                "at construction instead of passing steps=/simref= to run()."
            )
        super().initialize_run_parameters()

    # ------------------------------------------------------------------
    # Engine + sampling (one recompile per temperature stage)
    # ------------------------------------------------------------------
    def _make_engine(self) -> ThermalEngine:
        return ThermalEngine(
            self,
            rule=self.rule,
            upd_mode=self.upd_mode,
            T=float(self.T_schedule[0]),
            order=self.order,
            tie_flip_p=self.tie_flip_p,
        )

    def _sample_py(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        engine = self._make_engine()
        self._e_running = self.total_energy()
        n_stage = len(self.T_schedule)
        stage_e = np.empty(n_stage, dtype=np.float64)
        stage_m = np.empty(n_stage, dtype=np.float64)
        outer = (
            tqdm.tqdm(self.T_schedule, desc=type(self).__name__)
            if tqdm_on
            else self.T_schedule
        )
        for k, T_stage in enumerate(outer):
            self.T = float(T_stage)
            engine.T = float(T_stage)
            sweep = engine.compile_step()
            for _ in range(self.steps_per_T):
                self._record()
                self._e_running += sweep()
            stage_e[k] = self._e_running / float(self.N)
            stage_m[k] = self.magnetization()
        self.sa_temps = self.T_schedule.copy()
        self.sa_energy = stage_e
        self.sa_magn = stage_m

    # ------------------------------------------------------------------
    # Native pybind backend (sa kernel; guards in _pb_check_supported)
    # ------------------------------------------------------------------
    def _pb_check_supported(self) -> None:
        if self.field_mode == "spectral":
            raise NotImplementedError(
                "runlang='pb': the native kernel records the pairwise-only "
                "energy, so the spectral (TFCA) field would be dropped from "
                "the trace; field='spectral' runs on the python backend."
            )
        if self.schedule != "exponential":
            raise NotImplementedError(
                "runlang='pb': only schedule='exponential' matches the C "
                "cooling grid exactly (T_init·rate^k); the native "
                f"linear/logarithmic grids differ — schedule="
                f"{self.schedule!r} runs on the python backend."
            )
        if self.T_init <= 0.0:
            raise NotImplementedError(
                "runlang='pb': the native annealing loop stops at T <= 0, "
                "so an all-zero exponential schedule (T_init=0) cannot run; "
                "use runlang='py'."
            )
        self._pb_check_kinetics(
            rule=self.rule,
            order=self.order,
            tie_flip_p=self.tie_flip_p,
            zero_T=False,  # exponential stages are strictly positive
        )

    def _run_pb(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        """One native annealing run. ``T_final`` is passed as 0.0 to disable
        the C early-exit (``T > T_final`` loop guard): the python schedule
        generator runs ALL ``n_temperatures`` stages and ignores ``T_final``
        under 'exponential' (legacy semantics), so this reproduces the
        python-backend schedule exactly. Energies are the kernel's intensive
        pairwise values (== full Hamiltonian; the field is zero by guard)."""
        mod = load_ising_native()
        ni, nw, nptr = self._pb_csr()
        s_out, ene, magn, T_out = mod.sa_sampling(
            np.ascontiguousarray(self.s, dtype=np.int8),
            ni,
            nw,
            nptr,
            np.asarray(self.field, dtype=np.float64),
            float(self.T_init),
            0.0,  # no early exit — run the full python-equivalent schedule
            "exponential",
            float(self.cooling_rate),
            int(self.steps_per_T),
            int(self.n_temperatures),
            int(self.seed),
            PB_ORDER_MAP[self.order],
        )
        self.s = np.asarray(s_out, dtype=np.int8)
        ene = np.asarray(ene)
        magn = np.asarray(magn)
        self.ene = ene.tolist()
        self.magn = magn.tolist()
        self._e_running = self.total_energy()
        # End-of-stage curve: the trace is recorded BEFORE each sweep, so
        # stage k ends at the first record of stage k+1; the last stage
        # ends at the final state.
        k = self.steps_per_T
        self.sa_temps = np.asarray(T_out)
        self.sa_energy = np.append(ene[k::k], self.energy_intensive())
        self.sa_magn = np.append(magn[k::k], self.magnetization())

    def sa_curve(self) -> dict:
        """The per-stage annealing curve: end-of-stage intensive energy and
        magnetization at each schedule temperature."""
        if self.sa_energy is None:
            raise RuntimeError("No annealing run recorded yet; call run().")
        return {
            "T": self.sa_temps,
            "energy": self.sa_energy,
            "magn": self.sa_magn,
        }
