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

from .._thermal import THERMAL_RULES, TIE_FLIP_PRESETS, ThermalEngine
from .defaults import (
    ISING_MOVE_DEFAULT,
    ISING_MOVES,
    ISING_ORDER_DEFAULT,
    ISING_RULE_DEFAULT,
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
        Update schedule; Phase 1 wires ``'async'``.
    order : {'random', 'permutation', 'typewriter'}
        Site-visit order for async kinetics (default ``'random'`` — the
        legacy fixed-order loop survives as ``'typewriter'``).
    move : {'single', 'wolff', 'sw', 'spectral', 'exchange'}
        Elementary move; Phase 1 wires ``'single'``.
    tie_flip_p : float or str, optional
        Metropolis-only ΔE=0 acceptance: a probability in [0, 1] or a named
        preset (``'standard'`` 1.0 / ``'glauberT0'`` 0.5 / ``'frozen'`` 0.0).
        Rejected for ``glauber``/``heatbath`` (they bake their tie rate).
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
        if move != "single":
            raise NotImplementedError(
                f"move={move!r} is not wired on the new scheme classes yet "
                "(cluster/spectral/exchange moves land in later phases); use "
                "move='single' or the legacy IsingDynamics."
            )
        # Tie policy is a metropolis-only kinetic knob (plan §3b).
        if rule == "metropolis":
            if tie_flip_p is None:
                resolved_tie = ISING_TIE_FLIP_DEFAULT
            elif isinstance(tie_flip_p, str):
                try:
                    resolved_tie = TIE_FLIP_PRESETS[tie_flip_p]
                except KeyError:
                    raise ValueError(
                        f"Unknown tie_flip_p preset {tie_flip_p!r}; known "
                        f"presets: {tuple(TIE_FLIP_PRESETS)}."
                    ) from None
            else:
                resolved_tie = float(tie_flip_p)
        else:
            if tie_flip_p is not None:
                raise ValueError(
                    f"tie_flip_p is metropolis-only; rule {rule!r} bakes its "
                    "own ΔE=0 rate."
                )
            resolved_tie = None

        self.T = float(T)
        self.rule = rule
        self.upd_mode = upd_mode
        self.order = order
        self.move = move
        self.tie_flip_p = resolved_tie

        super().__init__(sg, **kwargs)

        # Setup-time capability check: constructing the engine validates
        # upd_mode/order/rule/T/tie combinations NOW, not mid-run.
        self._make_engine()

    def _make_engine(self) -> ThermalEngine:
        return ThermalEngine(
            self,
            rule=self.rule,
            upd_mode=self.upd_mode,
            T=self.T,
            order=self.order,
            tie_flip_p=self.tie_flip_p if self.rule == "metropolis" else None,
        )
