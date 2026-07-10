"""PottsMetropolis — fixed-temperature q-state Potts sampling (Layer-1 scheme).

The Potts scheme leaf on the shared spin-model engine (statsys unification
D-B8/D-B9): one system at one temperature, driven by the SAME
``statsys._thermal`` engines the Ising schemes use — no engine code is
rewritten, only the :class:`PottsBase` substrate differs. All Layer-2 choices
are parameters, never subclasses: acceptance ``rule``, update schedule
``upd_mode`` + async ``order``, elementary ``move``.

Capability boundaries (hard errors at setup, invariant #3 — never a silent
approximation):

- ``rule='heatbath'`` is q = 2 only: the shared kernel is the two-state
  conditional resample (M2/M4 coincidence); the true multi-state conditional
  is not wired. ``glauber`` stays valid for every q — the Barker acceptance
  ``1/(1 + e^{βΔE})`` satisfies detailed balance for any symmetric proposal.
- ``upd_mode='gillespie'`` is q = 2 only: the shared BKL engine (M3) assumes
  a deterministic local proposal; q > 2 needs per-(site,label) rate channels.
- Cluster moves (``wolff``/``sw``, M6/M5) require non-negative couplings: a
  joint relabel preserves LIKE pairs but not UNLIKE ones, so the satisfied
  state of an antiferromagnetic bond is not flip-covariant and the Ising
  detailed-balance bookkeeping does not carry over. The engines' ½-flip
  cluster resampling (each selected cluster gets a uniform random OTHER
  label) is a symmetric kernel on the FK conditional — uniform over cluster
  labels (M5) — so detailed balance holds; for q = 2 it is exactly M5's rule.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from .._thermal import (
    THERMAL_RULES,
    KawasakiEngine,
    SwendsenWangEngine,
    ThermalEngine,
    WolffEngine,
    resolve_tie_flip,
)
from ._native import load_potts_native
from .defaults import (
    POTTS_MOVE_DEFAULT,
    POTTS_MOVES,
    POTTS_OBS_SNAPSHOTS,
    POTTS_ORDER_DEFAULT,
    POTTS_Q_DEFAULT,
    POTTS_RULE_DEFAULT,
    POTTS_T_DEFAULT,
    POTTS_TIE_FLIP_DEFAULT,
    POTTS_UPD_MODE_DEFAULT,
)
from .PottsBase import PottsBase

if TYPE_CHECKING:
    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph

__all__ = ["PottsMetropolis"]

_TWO_STATE_ONLY_RULES = ("heatbath",)


class PottsMetropolis(PottsBase):
    """Fixed-T equilibrium Potts dynamics on the shared MCMC engine.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph carrying the couplings.
    q : int
        Number of Potts states (>= 2).
    T : float
        Temperature (``k_B = 1``); ``0`` runs the exact zero-temperature
        (quench) limit of the chosen rule.
    rule : {'metropolis', 'glauber', 'heatbath'}
        Acceptance rule (M1/M2/M4; ``heatbath`` is q = 2 only, see the
        module docstring).
    upd_mode : {'async', 'sync', 'gillespie'}
        Update schedule, ``move='single'`` only (``gillespie`` is q = 2
        only).
    order : {'random', 'permutation', 'typewriter'}
        Site-visit order for async kinetics (default ``'random'``; the
        legacy Python loop's fresh-permutation-per-sweep survives as
        ``'permutation'``). Applies to ``move='single'`` and
        ``move='exchange'`` under ``upd_mode='async'`` only.
    move : {'single', 'wolff', 'sw', 'exchange'}
        Elementary move: single-label change; Wolff single-cluster (M6) /
        Swendsen–Wang multi-cluster (M5) relabelings — non-negative
        couplings only, per-sweep cluster statistics in
        :attr:`cluster_sizes` / :attr:`cluster_counts`; ``'exchange'`` =
        Kawasaki pair swap (M8), conserves the per-state densities.
    tie_flip_p : float or str, optional
        Metropolis-only ΔE=0 acceptance: a probability in [0, 1] or a named
        preset (``'standard'`` 1.0 / ``'glauberT0'`` 0.5 / ``'frozen'``
        0.0). Rejected for ``glauber``/``heatbath`` and the cluster moves.
    **kwargs
        Forwarded to :class:`PottsBase` (``observables``, ``coupling_norm``,
        ``ic``, ``steps``/``simref``, ``seed``, ``savedisk``, ...).
    """

    def __init__(
        self,
        sg: "SignedGraph",
        q: int = POTTS_Q_DEFAULT,
        T: float = POTTS_T_DEFAULT,
        *,
        rule: str = POTTS_RULE_DEFAULT,
        upd_mode: str = POTTS_UPD_MODE_DEFAULT,
        order: str = POTTS_ORDER_DEFAULT,
        move: str = POTTS_MOVE_DEFAULT,
        tie_flip_p: float | str | None = None,
        **kwargs: Any,
    ) -> None:
        if rule not in THERMAL_RULES:
            raise ValueError(
                f"rule must be one of {THERMAL_RULES}; got {rule!r}."
            )
        if move not in POTTS_MOVES:
            raise ValueError(
                f"move must be one of {POTTS_MOVES}; got {move!r}."
            )
        if rule in _TWO_STATE_ONLY_RULES and int(q) > 2:
            raise NotImplementedError(
                f"rule={rule!r} with q={q}: the shared kernel is the "
                "two-state conditional resample (it coincides with 'glauber' "
                "only for q = 2); the multi-state heat-bath resample is not "
                "wired. Use rule='glauber' (valid for every q) or q = 2."
            )
        if upd_mode == "gillespie" and int(q) > 2:
            raise NotImplementedError(
                f"upd_mode='gillespie' with q={q}: the shared BKL engine "
                "(ref M3) assumes a deterministic local proposal, so it is "
                "two-state only; q > 2 needs per-(site,label) rate channels, "
                "which are not wired."
            )
        # Axis interactions (mirrors IsingMetropolis, plan §9): cluster moves
        # carry their own intrinsic schedule/acceptance, so the inapplicable
        # single-site axes are rejected — never silently ignored.
        if move in ("wolff", "sw"):
            if upd_mode != POTTS_UPD_MODE_DEFAULT:
                raise ValueError(
                    f"upd_mode does not apply to move={move!r} (the move "
                    "carries its own schedule); leave upd_mode at "
                    f"{POTTS_UPD_MODE_DEFAULT!r}."
                )
            if order != POTTS_ORDER_DEFAULT:
                raise ValueError(
                    f"order does not apply to move={move!r}; leave order at "
                    f"{POTTS_ORDER_DEFAULT!r}."
                )
            if rule != POTTS_RULE_DEFAULT:
                raise ValueError(
                    f"rule does not apply to move={move!r}: cluster moves "
                    "define their own (Fortuin-Kasteleyn) acceptance."
                )
            if tie_flip_p is not None:
                raise ValueError(f"tie_flip_p does not apply to move={move!r}.")
            resolved_tie = None
        else:
            if move == "exchange" and upd_mode != POTTS_UPD_MODE_DEFAULT:
                raise ValueError(
                    "move='exchange' has asynchronous kinetics only (the "
                    "swap-attempt sequence is its schedule); leave upd_mode "
                    f"at {POTTS_UPD_MODE_DEFAULT!r}."
                )
            # Tie policy is a metropolis-only kinetic knob (plan §3b); the
            # resolution logic is shared by all thermal schemes (_thermal).
            resolved_tie = resolve_tie_flip(
                rule, tie_flip_p, default=POTTS_TIE_FLIP_DEFAULT
            )

        self.rule = rule
        self.upd_mode = upd_mode
        self.order = order
        self.move = move
        self.tie_flip_p = resolved_tie
        self.cluster_sizes: list[int] | None = None
        self.cluster_counts: list[int] | None = None

        super().__init__(sg, q=q, T=T, **kwargs)

        # Setup-time capability check: constructing the engine validates
        # upd_mode/order/rule/T/tie combinations NOW, not mid-run. The
        # cluster engines additionally need the couplings (sign check), so
        # for those moves the same validation happens in _make_engine on
        # first use, still before any sweep runs.
        if move == "single":
            self._make_engine()

    # ------------------------------------------------------------------
    # Run-dirname schema (Phase C)
    # ------------------------------------------------------------------
    def _physics_tokens(self) -> list:
        return super()._physics_tokens()

    def _axis_tokens(self) -> list:
        cluster = self.move in ("wolff", "sw")
        return [
            ("rule", None if cluster else self.rule, POTTS_RULE_DEFAULT),
            ("move", self.move, POTTS_MOVE_DEFAULT),
            ("upd", self.upd_mode, POTTS_UPD_MODE_DEFAULT),
            (
                "ord",
                self.order if self.upd_mode == "async" else None,
                POTTS_ORDER_DEFAULT,
            ),
            (
                "tf",
                (
                    None
                    if cluster or self.rule != "metropolis"
                    else self.tie_flip_p
                ),
                1.0,
            ),
        ]

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
        if self.move in ("wolff", "sw"):
            if np.any(self._nbr_J < 0.0):
                raise NotImplementedError(
                    f"move={self.move!r} needs non-negative couplings: a "
                    "joint relabel does not preserve an unlike pair, so "
                    "antiferromagnetic Potts bonds break the FK "
                    "detailed-balance bookkeeping (unlike the Ising joint "
                    "flip). Use a single-site move on signed graphs."
                )
            if self.move == "wolff":
                return WolffEngine(self, T=self.T)
            return SwendsenWangEngine(self, T=self.T)
        # exchange
        return KawasakiEngine(
            self,
            rule=self.rule,
            T=self.T,
            order=self.order,
            tie_flip_p=(self.tie_flip_p if self.rule == "metropolis" else None),
        )

    # ------------------------------------------------------------------
    # Vectorized sync backend (np/cu) — capability gate + engine factory
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # Native pybind backend (potts_sampling kernel; guards below)
    # ------------------------------------------------------------------
    def _pb_check_supported(self) -> None:
        """What the compiled kernel (``LRGSG_potts.c``) cannot represent is a
        hard capability error (invariant #3): it implements the Metropolis
        acceptance only, visits sites uniformly at random WITH replacement
        (no other order), and accepts every ΔE <= 0 proposal — so its tie
        policy is 1.0 at EVERY temperature, including T = 0."""
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
        if self.rule != "metropolis":
            raise NotImplementedError(
                "runlang='pb': the native kernel implements the metropolis "
                f"acceptance only; rule={self.rule!r} is python-only."
            )
        if self.order != "random":
            raise NotImplementedError(
                "runlang='pb': the native kernel visits sites uniformly at "
                f"random (with replacement) only; order={self.order!r} is "
                "python-only."
            )
        if self.tie_flip_p is not None and float(self.tie_flip_p) != 1.0:
            raise NotImplementedError(
                "runlang='pb': the native kernel accepts every ΔE <= 0 "
                "proposal (tie policy 1.0 at every temperature, including "
                f"T = 0); tie_flip_p={self.tie_flip_p} is python-only."
            )
        if POTTS_OBS_SNAPSHOTS in self._selected_obs:
            raise NotImplementedError(
                "runlang='pb': the native kernel returns energy/magn traces "
                "and the final state only; deselect 'snapshots' or use "
                "runlang='py'."
            )

    def _run_pb(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        """One fixed-T native run: same J (coupling_norm applied), intensive
        energies. The kernel records AFTER each sweep, so the §3b
        record-then-sweep trace is reconstructed exactly by prepending the
        initial values and dropping the kernel's last record."""
        mod = load_potts_native()
        ni, nw, nptr = self._pb_csr()
        e0 = self.energy_intensive()
        m0 = self.order_parameter()
        s_out, ene, magn = mod.potts_sampling(
            np.ascontiguousarray(self.s, dtype=np.int32),
            ni,
            nw,
            nptr,
            int(self.q),
            float(self.T),
            int(self.steps),
            int(self.seed),
            True,
        )
        self.s = np.asarray(s_out, dtype=np.int32)
        self.ene = [e0] + (
            np.asarray(ene, dtype=np.float64)[:-1] / float(self.N)
        ).tolist()
        self.magn = [m0] + np.asarray(magn, dtype=np.float64)[:-1].tolist()
        self._e_running = self.total_energy()

    def _sample_py(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        super()._sample_py(tqdm_on=tqdm_on, verbose=verbose)
        # Per-move extras, mirrored off the engine after the run.
        if self.move in ("wolff", "sw"):
            self.cluster_sizes = getattr(self._engine, "cluster_sizes", None)
            self.cluster_counts = getattr(self._engine, "cluster_counts", None)
