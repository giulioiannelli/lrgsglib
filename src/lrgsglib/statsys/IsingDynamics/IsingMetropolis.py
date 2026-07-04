"""IsingMetropolis — fixed-temperature equilibrium sampling (Layer-1 scheme).

The first thin scheme leaf of the Ising refactor (plan §2/§4): one system at
one temperature, single-site kinetics driven by the shared
``statsys._thermal.ThermalEngine``. All Layer-2 choices are parameters, never
subclasses: acceptance ``rule`` (metropolis M1 / glauber M2 / heatbath M4 —
citations in ``.agents/ref/00-REFERENCES.md`` §MCMC), update schedule
``upd_mode`` + async ``order``, elementary ``move``.

Phase-1 slice: ``move='single'`` and ``upd_mode='async'`` on the Python
backend. The other axis values are validated here at construction and raise a
hard capability error (invariant #3: never a silent fallback); they are wired
in later phases (cluster moves, sync/gillespie, native backends).
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
        Update schedule (single-site moves only).
    order : {'random', 'permutation', 'typewriter'}
        Site-visit order for async kinetics (default ``'random'`` — the
        legacy fixed-order loop survives as ``'typewriter'``). Applies to
        ``move='single'`` and ``move='exchange'``.
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
