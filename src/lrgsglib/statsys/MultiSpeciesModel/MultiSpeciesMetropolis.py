"""MultiSpeciesMetropolis — fixed-temperature k-species sampling (Layer-1).

The multi-species scheme leaf on the shared spin-model engine (statsys
unification D-B8/D-B9, closing the spin quad): one system at one
temperature, driven by the SAME ``statsys._thermal`` engines the other spin
schemes use — no engine code is rewritten, only the
:class:`MultiSpeciesBase` substrate differs (a (component, label) proposal
on a k-column integer state). All Layer-2 choices are parameters, never
subclasses: acceptance ``rule``, update schedule ``upd_mode`` + async
``order``, elementary ``move``.

Capability boundaries (hard errors at setup, invariant #3 — never a silent
approximation):

- ``rule='heatbath'`` is not wired: the conditional resample over the
  coupled (component, label) channels is not implemented. ``glauber``
  (Barker) stays valid — the (component, label) proposal is symmetric.
- ``upd_mode='gillespie'`` is not wired: the shared BKL engine (M3) assumes
  a deterministic local proposal; the multi-species proposal draws a random
  component AND label (per-channel rate tables are not wired).
- Cluster moves are not wired: the FK construction for the COUPLED
  multi-species Hamiltonian (inter-species delta terms) does not reduce to
  the shared single-label bookkeeping.
- The species-interaction matrix must be SYMMETRIC (validated by the base).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

import numpy as np

from .._thermal import (
    THERMAL_RULES,
    KawasakiEngine,
    ThermalEngine,
    resolve_tie_flip,
)
from ._native import load_multispec_native
from .defaults import (
    MULTISPEC_MOVE_DEFAULT,
    MULTISPEC_MOVES,
    MULTISPEC_OBS_MAGN,
    MULTISPEC_OBS_SNAPSHOTS,
    MULTISPEC_ORDER_DEFAULT,
    MULTISPEC_Q_DEFAULT,
    MULTISPEC_RULE_DEFAULT,
    MULTISPEC_SPECIES_DEFAULT,
    MULTISPEC_T_DEFAULT,
    MULTISPEC_TIE_FLIP_DEFAULT,
    MULTISPEC_UPD_MODE_DEFAULT,
)
from .MultiSpeciesBase import MultiSpeciesBase

if TYPE_CHECKING:
    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph

__all__ = ["MultiSpeciesMetropolis"]

_UNWIRED_RULES = ("heatbath",)
_UNWIRED_MOVES = ("wolff", "sw")


class MultiSpeciesMetropolis(MultiSpeciesBase):
    """Fixed-T equilibrium multi-species dynamics on the shared MCMC engine.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph carrying the couplings.
    species : int
        Number of coupled Potts components per site (k >= 1).
    q_per_species : int or sequence of int
        Number of states per component (each >= 2).
    T : float
        Temperature (``k_B = 1``); ``0`` runs the exact zero-temperature
        (quench) limit of the chosen rule.
    rule : {'metropolis', 'glauber'}
        Acceptance rule (M1/M2; ``heatbath`` is not wired, see the module
        docstring).
    upd_mode : {'async', 'sync'}
        Update schedule, ``move='single'`` only (``gillespie`` is not
        wired).
    order : {'random', 'permutation', 'typewriter'}
        Site-visit order for async kinetics (default ``'random'``; the
        legacy Python loop's fresh-permutation-per-sweep survives as
        ``'permutation'``). Applies to ``move='single'`` and
        ``move='exchange'`` under ``upd_mode='async'`` only.
    move : {'single', 'exchange'}
        Elementary move: single-component relabel, or ``'exchange'`` =
        Kawasaki FULL-site swap (M8), conserving every per-species density.
        Cluster moves are not wired for the coupled Hamiltonian.
    tie_flip_p : float or str, optional
        Metropolis-only ΔE=0 acceptance: a probability in [0, 1] or a named
        preset (``'standard'`` 1.0 / ``'glauberT0'`` 0.5 / ``'frozen'``
        0.0). Rejected for ``glauber``.
    **kwargs
        Forwarded to :class:`MultiSpeciesBase` (``interaction_matrix``,
        ``observables``, ``coupling_norm``, ``ic``, ``steps``/``simref``,
        ``seed``, ``savedisk``, ...).
    """

    def __init__(
        self,
        sg: "SignedGraph",
        species: int = MULTISPEC_SPECIES_DEFAULT,
        q_per_species: int | Any = MULTISPEC_Q_DEFAULT,
        T: float = MULTISPEC_T_DEFAULT,
        *,
        rule: str = MULTISPEC_RULE_DEFAULT,
        upd_mode: str = MULTISPEC_UPD_MODE_DEFAULT,
        order: str = MULTISPEC_ORDER_DEFAULT,
        move: str = MULTISPEC_MOVE_DEFAULT,
        tie_flip_p: float | str | None = None,
        **kwargs: Any,
    ) -> None:
        if rule not in THERMAL_RULES:
            raise ValueError(
                f"rule must be one of {THERMAL_RULES}; got {rule!r}."
            )
        if move in _UNWIRED_MOVES:
            raise NotImplementedError(
                f"move={move!r}: cluster moves are not wired for the coupled "
                "multi-species Hamiltonian — the inter-species delta terms "
                "do not reduce to the shared single-label FK bookkeeping. "
                f"Available moves: {MULTISPEC_MOVES}."
            )
        if move not in MULTISPEC_MOVES:
            raise ValueError(
                f"move must be one of {MULTISPEC_MOVES}; got {move!r}."
            )
        if rule in _UNWIRED_RULES:
            raise NotImplementedError(
                f"rule={rule!r}: the conditional resample over the coupled "
                "(component, label) channels is not wired. Use "
                "rule='glauber' (valid for the symmetric proposal) or "
                "rule='metropolis'."
            )
        if upd_mode == "gillespie":
            raise NotImplementedError(
                "upd_mode='gillespie': the shared BKL engine (ref M3) "
                "assumes a deterministic local proposal; the multi-species "
                "proposal draws a random component AND label (per-channel "
                "rate tables are not wired)."
            )
        if move == "exchange" and upd_mode != MULTISPEC_UPD_MODE_DEFAULT:
            raise ValueError(
                "move='exchange' has asynchronous kinetics only (the "
                "swap-attempt sequence is its schedule); leave upd_mode "
                f"at {MULTISPEC_UPD_MODE_DEFAULT!r}."
            )
        # Tie policy is a metropolis-only kinetic knob (plan §3b); the
        # resolution logic is shared by all thermal schemes (_thermal).
        resolved_tie = resolve_tie_flip(
            rule, tie_flip_p, default=MULTISPEC_TIE_FLIP_DEFAULT
        )

        self.rule = rule
        self.upd_mode = upd_mode
        self.order = order
        self.move = move
        self.tie_flip_p = resolved_tie

        super().__init__(
            sg, species=species, q_per_species=q_per_species, T=T, **kwargs
        )

        # Setup-time capability check: constructing the engine validates
        # upd_mode/order/rule/T/tie combinations NOW, not mid-run.
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
            ("rule", None if cluster else self.rule, MULTISPEC_RULE_DEFAULT),
            ("move", self.move, MULTISPEC_MOVE_DEFAULT),
            ("upd", self.upd_mode, MULTISPEC_UPD_MODE_DEFAULT),
            (
                "ord",
                self.order if self.upd_mode == "async" else None,
                MULTISPEC_ORDER_DEFAULT,
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
        # exchange
        self._ensure_couplings()
        return KawasakiEngine(
            self,
            rule=self.rule,
            T=self.T,
            order=self.order,
            tie_flip_p=(self.tie_flip_p if self.rule == "metropolis" else None),
        )

    # ------------------------------------------------------------------
    # Native pybind backend (multispec_sampling kernel; guards below)
    # ------------------------------------------------------------------
    def _pb_check_supported(self) -> None:
        """What the compiled kernel (``LRGSG_multispec.c``) cannot represent
        is a hard capability error (invariant #3): it implements the
        Metropolis acceptance only, visits sites uniformly at random WITH
        replacement, accepts every ΔE <= 0 proposal (tie policy 1.0 at
        every temperature, including T = 0), assumes the IDENTITY
        interaction matrix and a UNIFORM q, and records the energy only."""
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
        if not np.array_equal(self.interaction_matrix, np.eye(self.species)):
            raise NotImplementedError(
                "runlang='pb': the native kernel assumes the IDENTITY "
                "species-interaction matrix; a custom interaction_matrix "
                "is python-only."
            )
        if len(set(self._q_per_species)) != 1:
            raise NotImplementedError(
                "runlang='pb': the native kernel assumes a UNIFORM number "
                f"of states per species; q_per_species={self._q_per_species} "
                "is python-only."
            )
        if MULTISPEC_OBS_MAGN in self._selected_obs:
            raise NotImplementedError(
                "runlang='pb': the native kernel records the energy only; "
                "deselect 'magn' (observables=('energy',)) or use "
                "runlang='py'."
            )
        if MULTISPEC_OBS_SNAPSHOTS in self._selected_obs:
            raise NotImplementedError(
                "runlang='pb': the native kernel returns the energy trace "
                "and the final state only; deselect 'snapshots' or use "
                "runlang='py'."
            )

    def _run_pb(self, tqdm_on: bool = False, verbose: bool = False) -> None:
        """One fixed-T native run: same J (coupling_norm applied), intensive
        energies. The kernel records AFTER each sweep, so the §3b
        record-then-sweep trace is reconstructed exactly by prepending the
        initial value and dropping the kernel's last record."""
        mod = load_multispec_native()
        ni, nw, nptr = self._pb_csr()
        e0 = self.energy_intensive()
        s_out, ene = mod.multispec_sampling(
            np.ascontiguousarray(self.s, dtype=np.int32).reshape(-1),
            ni,
            nw,
            nptr,
            int(self.N),
            int(self.species),
            int(self._q_per_species[0]),
            float(self.T),
            int(self.steps),
            int(self.seed),
            True,
        )
        self.s = np.asarray(s_out, dtype=np.int32).reshape(self.N, self.species)
        self.ene = [e0] + (
            np.asarray(ene, dtype=np.float64)[:-1] / float(self.N)
        ).tolist()
        self._e_running = self.total_energy()
