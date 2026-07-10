"""IsingMetropolis — fixed-temperature equilibrium sampling (Layer-1 scheme).

The first thin scheme leaf of the Ising refactor (plan §2/§4): one system at
one temperature, single-site kinetics driven by the shared
``statsys._thermal.ThermalEngine``. All Layer-2 choices are parameters, never
subclasses: acceptance ``rule`` (metropolis M1 / glauber M2 / heatbath M4 —
citations in ``.agents/ref/00-REFERENCES.md`` §MCMC), update schedule
``upd_mode`` + async ``order``, elementary ``move``.

Phase-2 state: every move (single/wolff/sw/spectral/exchange) and every
update mode (async/sync/gillespie, single-site moves only) is wired on the
Python backend. Inapplicable axis combinations are validated here at
construction and raise a hard capability error (invariant #3: never a silent
fallback); the native backends land next.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from .._spectral import (
    SpectralMetropolisEngine,
    best_eigenvector_seed,
    eigvec_subspace_matrix,
    rbim_subspace_energies,
)
from .._thermal import (
    THERMAL_RULES,
    KawasakiEngine,
    SwendsenWangEngine,
    ThermalEngine,
    WolffEngine,
    resolve_tie_flip,
)
from ._native import PB_ORDER_MAP, load_ising_native
from .defaults import (
    ISING_MOVE_DEFAULT,
    ISING_MOVES,
    ISING_ORDER_DEFAULT,
    ISING_RULE_DEFAULT,
    ISING_SPECTRAL_CHUNK_SIZE_DEFAULT,
    ISING_SPECTRAL_N_MODES_DEFAULT,
    ISING_SPECTRAL_POLISH_DEFAULT,
    ISING_SPECTRAL_POLISH_SWEEPS_DEFAULT,
    ISING_SPECTRAL_SIGMA_INIT_DEFAULT,
    ISING_T_DEFAULT,
    ISING_TIE_FLIP_DEFAULT,
    ISING_UPD_MODE_DEFAULT,
)
from .IsingBase import IsingBase

if TYPE_CHECKING:
    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph

__all__ = ["IsingMetropolis"]


class IsingMetropolis(IsingBase):
    """Fixed-T equilibrium Ising dynamics on the shared MCMC engine.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph carrying the couplings.
    T : float
        Temperature (``k_B = 1``); ``0`` runs the exact zero-temperature
        (quench) limit of the chosen rule.
    rule : {'metropolis', 'glauber', 'heatbath'}
        Acceptance rule (module docstring for citations).
    upd_mode : {'async', 'sync', 'gillespie'}
        Update schedule, ``move='single'`` only: random-sequential attempts;
        sublattice-parallel synchronous updates (greedy graph coloring, M4);
        rejection-free continuous-time BKL kinetics (M3, one recorded step =
        one unit of time = one sweep-equivalent).
    order : {'random', 'permutation', 'typewriter'}
        Site-visit order for async kinetics (default ``'random'`` — the
        legacy fixed-order loop survives as ``'typewriter'``). Applies to
        ``move='single'`` and ``move='exchange'`` under ``upd_mode='async'``
        only.
    move : {'single', 'wolff', 'sw', 'spectral', 'exchange'}
        Elementary move: single-spin flip; Wolff single-cluster (M6) /
        Swendsen–Wang multi-cluster (M5) — zero-field only, per-sweep
        cluster statistics in :attr:`cluster_sizes` / :attr:`cluster_counts`;
        ``'spectral'`` = Metropolis in the eigenvector-coefficient space
        (extras in :attr:`spectral_coeffs` / :attr:`spectral_best_spins` /
        :attr:`spectral_best_energy`, intensive); ``'exchange'`` = Kawasaki
        pair swap (M8), conserves the magnetization.
    tie_flip_p : float or str, optional
        Metropolis-only ΔE=0 acceptance: a probability in [0, 1] or a named
        preset (``'standard'`` 1.0 / ``'glauberT0'`` 0.5 / ``'frozen'`` 0.0).
        Rejected for ``glauber``/``heatbath`` (they bake their tie rate) and
        for the cluster/spectral moves (the axis does not apply).
    spectral_n_modes, spectral_sigma_init, spectral_chunk_size,
    spectral_polish, spectral_polish_sweeps
        Parameters of the spectral move (``move='spectral'`` only).
    **kwargs
        Forwarded to :class:`IsingBase` (``observables``, ``coupling_norm``,
        ``field``, ``ic``, ``steps``/``simref``, ``seed``, ``savedisk``, ...).
    """

    def __init__(
        self,
        sg: "SignedGraph",
        T: float = ISING_T_DEFAULT,
        *,
        rule: str = ISING_RULE_DEFAULT,
        upd_mode: str = ISING_UPD_MODE_DEFAULT,
        order: str = ISING_ORDER_DEFAULT,
        move: str = ISING_MOVE_DEFAULT,
        tie_flip_p: float | str | None = None,
        spectral_n_modes: int = ISING_SPECTRAL_N_MODES_DEFAULT,
        spectral_sigma_init: float = ISING_SPECTRAL_SIGMA_INIT_DEFAULT,
        spectral_chunk_size: int = ISING_SPECTRAL_CHUNK_SIZE_DEFAULT,
        spectral_polish: bool = ISING_SPECTRAL_POLISH_DEFAULT,
        spectral_polish_sweeps: int = ISING_SPECTRAL_POLISH_SWEEPS_DEFAULT,
        **kwargs: Any,
    ) -> None:
        if rule not in THERMAL_RULES:
            raise ValueError(
                f"rule must be one of {THERMAL_RULES}; got {rule!r}."
            )
        if move not in ISING_MOVES:
            raise ValueError(
                f"move must be one of {ISING_MOVES}; got {move!r}."
            )
        # Axis interactions (plan §9): cluster and spectral moves carry
        # their own intrinsic schedule/acceptance, so the inapplicable
        # single-site axes are rejected — never silently ignored.
        if move in ("wolff", "sw", "spectral"):
            if upd_mode != ISING_UPD_MODE_DEFAULT:
                raise ValueError(
                    f"upd_mode does not apply to move={move!r} (the move "
                    "carries its own schedule); leave upd_mode at "
                    f"{ISING_UPD_MODE_DEFAULT!r}."
                )
            if order != ISING_ORDER_DEFAULT:
                raise ValueError(
                    f"order does not apply to move={move!r}; leave order at "
                    f"{ISING_ORDER_DEFAULT!r}."
                )
        if move == "exchange" and upd_mode != ISING_UPD_MODE_DEFAULT:
            raise ValueError(
                "move='exchange' has asynchronous kinetics only (the swap-"
                "attempt sequence is its schedule); leave upd_mode at "
                f"{ISING_UPD_MODE_DEFAULT!r}."
            )
        if move in ("wolff", "sw"):
            if rule != ISING_RULE_DEFAULT:
                raise ValueError(
                    f"rule does not apply to move={move!r}: cluster moves "
                    "define their own (Fortuin-Kasteleyn) acceptance."
                )
            if tie_flip_p is not None:
                raise ValueError(f"tie_flip_p does not apply to move={move!r}.")
            resolved_tie = None
        elif move == "spectral":
            if rule != "metropolis":
                raise ValueError(
                    "move='spectral' is a Metropolis walk in coefficient "
                    f"space; rule {rule!r} does not apply."
                )
            if tie_flip_p is not None:
                raise ValueError(
                    "tie_flip_p does not apply to move='spectral' "
                    "(continuous proposals; ties have measure zero)."
                )
            resolved_tie = None
        else:
            # Tie policy is a metropolis-only kinetic knob (plan §3b); the
            # resolution logic is shared by all thermal schemes (_thermal).
            resolved_tie = resolve_tie_flip(
                rule, tie_flip_p, default=ISING_TIE_FLIP_DEFAULT
            )

        self.T = float(T)
        self.rule = rule
        self.upd_mode = upd_mode
        self.order = order
        self.move = move
        self.tie_flip_p = resolved_tie
        self.spectral_n_modes = int(spectral_n_modes)
        self.spectral_sigma_init = float(spectral_sigma_init)
        self.spectral_chunk_size = int(spectral_chunk_size)
        self.spectral_polish = bool(spectral_polish)
        self.spectral_polish_sweeps = int(spectral_polish_sweeps)
        self.cluster_sizes: list[int] | None = None
        self.cluster_counts: list[int] | None = None
        self.spectral_coeffs: np.ndarray | None = None
        self.spectral_best_spins: np.ndarray | None = None
        self.spectral_best_energy: float | None = None

        super().__init__(sg, **kwargs)

        # Setup-time capability check: constructing the engine validates
        # upd_mode/order/rule/T/tie combinations NOW, not mid-run. The
        # cluster/spectral engines additionally need the couplings (and the
        # spectral subspace), which exist only after init — for those moves
        # the same validation happens in _make_engine on first use, still
        # before any sweep runs.
        if move == "single":
            self._make_engine()
        elif move in ("wolff", "sw") and np.any(
            np.asarray(self.field, dtype=np.float64) != 0.0
        ):
            raise NotImplementedError(
                f"move={move!r} with a nonzero external field needs the "
                "ghost-spin extension, which is not wired; use h = 0 or a "
                "single-site move."
            )

    def _make_engine(self):
        if self.move == "single":
            return ThermalEngine(
                self,
                rule=self.rule,
                upd_mode=self.upd_mode,
                T=self.T,
                order=self.order,
                tie_flip_p=(
                    self.tie_flip_p if self.rule == "metropolis" else None
                ),
            )
        self._ensure_couplings()
        if self.move == "wolff":
            return WolffEngine(self, T=self.T)
        if self.move == "sw":
            return SwendsenWangEngine(self, T=self.T)
        if self.move == "exchange":
            return KawasakiEngine(
                self,
                rule=self.rule,
                T=self.T,
                order=self.order,
                tie_flip_p=(
                    self.tie_flip_p if self.rule == "metropolis" else None
                ),
            )
        # spectral
        n_modes = min(self.spectral_n_modes, self.N)
        V = eigvec_subspace_matrix(self.sg, n_modes)
        rbim_E = rbim_subspace_energies(self.sg, n_modes)
        seed_spins, seed_mode = best_eigenvector_seed(self.sg, n_modes)
        return SpectralMetropolisEngine(
            self,
            T=self.T,
            subspace_vecs=V,
            rbim_energies=rbim_E,
            seed_spins=seed_spins,
            seed_mode=seed_mode,
            sigma_init=self.spectral_sigma_init,
            chunk_size=self.spectral_chunk_size,
            polish=self.spectral_polish,
            polish_sweeps=self.spectral_polish_sweeps,
            nbr_idx=self._nbr_idx,
            nbr_J=self._nbr_J,
            nbr_ptr=self._nbr_ptr,
            field=np.asarray(self.field, dtype=np.float64),
        )

    # ------------------------------------------------------------------
    # Native pybind backend (met kernel; guards in _pb_check_supported)
    # ------------------------------------------------------------------
    def _pb_check_supported(self) -> None:
        if self.move != "single":
            raise NotImplementedError(
                "runlang='pb' is wired for move='single' only; "
                f"move={self.move!r} runs on the python backend."
            )
        if self.upd_mode != "async":
            raise NotImplementedError(
                "runlang='pb': the native kernel is asynchronous only; "
                f"upd_mode={self.upd_mode!r} runs on the python backend."
            )
        self._pb_check_kinetics(
            rule=self.rule,
            order=self.order,
            tie_flip_p=self.tie_flip_p,
            zero_T=(self.T == 0.0),
        )

    def _run_pb(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        """One fixed-T native run: same J (coupling_norm applied), same
        record-then-sweep cadence, intensive energies (the kernel's
        ``calc_totEnergy`` is already E/N; field is zero by the guard, so
        pairwise == full Hamiltonian)."""
        mod = load_ising_native()
        ni, nw, nptr = self._pb_csr()
        s_out, ene, magn = mod.metropolis_sampling(
            np.ascontiguousarray(self.s, dtype=np.int8),
            ni,
            nw,
            nptr,
            np.asarray(self.field, dtype=np.float64),
            float(self.T),
            int(self.steps),
            int(self.seed),
            PB_ORDER_MAP[self.order],
        )
        self.s = np.asarray(s_out, dtype=np.int8)
        self.ene = np.asarray(ene).tolist()
        self.magn = np.asarray(magn).tolist()
        self._e_running = self.total_energy()

    # ------------------------------------------------------------------
    # Run-dirname schema (Phase C)
    # ------------------------------------------------------------------
    def _physics_tokens(self) -> list:
        return [("T", float(self.T)), self._field_token()]

    def _axis_tokens(self) -> list:
        cluster = self.move in ("wolff", "sw")
        # tie-flip policy default: accept ties at T>0, reject at T=0
        tf_default = 1.0 if self.T > 0.0 else 0.0
        return [
            ("rule", None if cluster else self.rule, ISING_RULE_DEFAULT),
            ("move", self.move, ISING_MOVE_DEFAULT),
            ("upd", self.upd_mode, ISING_UPD_MODE_DEFAULT),
            (
                "ord",
                self.order if self.upd_mode == "async" else None,
                ISING_ORDER_DEFAULT,
            ),
            (
                "tf",
                (
                    None
                    if cluster or self.rule != "metropolis"
                    else self.tie_flip_p
                ),
                tf_default,
            ),
        ]

    def _np_check_supported(self) -> None:
        """The vectorized backends implement the SYNC schedule only —
        true async is inherently sequential (python/pb keep it)."""
        if self.upd_mode != "sync":
            raise NotImplementedError(
                "np/cu backends vectorize upd_mode='sync' (sublattice-"
                f"parallel); upd_mode={self.upd_mode!r} is sequential — "
                "use runlang='py' or 'pb'."
            )
        if self.move != "single":
            raise NotImplementedError(
                f"np/cu backends implement move='single'; move="
                f"{self.move!r} is python-only."
            )
        if self.rule not in ("metropolis", "glauber"):
            raise NotImplementedError(
                f"np/cu backends implement metropolis/glauber; rule="
                f"{self.rule!r} is python-only."
            )

    def _make_vec_engine(self, xp):
        from .._vectorized import VectorSyncEngine

        return VectorSyncEngine(
            self,
            rule=self.rule,
            T=self.T,
            tie_flip_p=(self.tie_flip_p if self.rule == "metropolis" else None),
            xp=xp,
        )

    def _sample_py(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        super()._sample_py(tqdm_on=tqdm_on, verbose=verbose)
        # Per-move extras, mirrored off the engine after the run.
        engine = self._engine
        if self.move in ("wolff", "sw"):
            self.cluster_sizes = getattr(engine, "cluster_sizes", None)
            self.cluster_counts = getattr(engine, "cluster_counts", None)
        elif self.move == "spectral":
            self.spectral_coeffs = engine.coeffs
            self.spectral_best_spins = engine.best_spins
            self.spectral_best_energy = engine.best_energy / float(self.N)
