"""XYMetropolis — fixed-temperature planar-rotator sampling (Layer-1 scheme).

The XY scheme leaf on the shared spin-model engine (statsys unification
D-B8/D-B9): one system at one temperature, driven by the SAME
``statsys._thermal`` engines the Ising and Potts schemes use — no engine code
is rewritten, only the :class:`XYBase` substrate differs (the first
continuous-state substrate). All Layer-2 choices are parameters, never
subclasses: acceptance ``rule``, update schedule ``upd_mode`` + async
``order``, elementary ``move``, proposal width ``delta``.

Capability boundaries (hard errors at setup, invariant #3 — never a silent
approximation):

- ``rule='heatbath'`` is not wired: for a continuous state the heat-bath
  update is a conditional resample from the local von-Mises distribution,
  which is not implemented. ``glauber`` (Barker) stays valid — the
  acceptance ``1/(1 + e^{βΔE})`` satisfies detailed balance for any
  symmetric proposal, and the ±delta perturbation is symmetric.
- ``upd_mode='gillespie'`` is not wired: the shared BKL engine (M3) assumes
  a deterministic local proposal; a continuous proposal is a fresh random
  draw, so single-proposal rates do not define the kinetics (an angular
  rate DENSITY would be needed).
- ``tie_flip_p`` is rejected: ``ΔE = 0`` has measure zero for a continuous
  state, so the tie knob has no XY meaning.
- Cluster moves (``wolff``/``sw``) use M6's embedded-Ising reflection and
  stay valid on SIGNED couplings (a reflection is an involution, so the
  embedded bond satisfaction is flip-covariant exactly as in Ising —
  contrast Potts, where the joint relabel forces a non-negative guard).
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
)
from ._native import load_xy_native
from .defaults import (
    XY_DELTA_DEFAULT,
    XY_MOVE_DEFAULT,
    XY_MOVES,
    XY_OBS_SNAPSHOTS,
    XY_ORDER_DEFAULT,
    XY_RULE_DEFAULT,
    XY_T_DEFAULT,
    XY_UPD_MODE_DEFAULT,
)
from .XYBase import XYBase

if TYPE_CHECKING:
    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph

__all__ = ["XYMetropolis"]

_DISCRETE_ONLY_RULES = ("heatbath",)


class XYMetropolis(XYBase):
    """Fixed-T equilibrium XY dynamics on the shared MCMC engine.

    Parameters
    ----------
    sg : SignedGraph
        The signed graph carrying the couplings.
    T : float
        Temperature (``k_B = 1``); ``0`` runs the exact zero-temperature
        (quench) limit of the chosen rule.
    delta : float, optional
        Maximum angular perturbation of the single-site proposal (radians);
        applies to ``move='single'`` only (the cluster and exchange moves
        have no proposal width) and defaults to the historical 0.5.
    rule : {'metropolis', 'glauber'}
        Acceptance rule (M1/M2; ``heatbath`` is not wired for a continuous
        state, see the module docstring).
    upd_mode : {'async', 'sync'}
        Update schedule, ``move='single'`` only (``gillespie`` is not
        available for a continuous state).
    order : {'random', 'permutation', 'typewriter'}
        Site-visit order for async kinetics (default ``'random'``; the
        legacy Python loop's fresh-permutation-per-sweep survives as
        ``'permutation'``). Applies to ``move='single'`` and
        ``move='exchange'`` under ``upd_mode='async'`` only.
    move : {'single', 'wolff', 'sw', 'exchange'}
        Elementary move: single-angle perturbation; Wolff single-cluster /
        Swendsen–Wang multi-cluster reflections via M6's embedded-Ising
        construction — valid on signed couplings, per-sweep cluster
        statistics in :attr:`cluster_sizes` / :attr:`cluster_counts`;
        ``'exchange'`` = Kawasaki pair swap (M8), conserves the multiset of
        angles.
    **kwargs
        Forwarded to :class:`XYBase` (``observables``, ``coupling_norm``,
        ``ic``, ``steps``/``simref``, ``seed``, ``savedisk``, ...).
    """

    def __init__(
        self,
        sg: "SignedGraph",
        T: float = XY_T_DEFAULT,
        delta: float | None = None,
        *,
        rule: str = XY_RULE_DEFAULT,
        upd_mode: str = XY_UPD_MODE_DEFAULT,
        order: str = XY_ORDER_DEFAULT,
        move: str = XY_MOVE_DEFAULT,
        tie_flip_p: float | str | None = None,
        **kwargs: Any,
    ) -> None:
        if rule not in THERMAL_RULES:
            raise ValueError(
                f"rule must be one of {THERMAL_RULES}; got {rule!r}."
            )
        if move not in XY_MOVES:
            raise ValueError(f"move must be one of {XY_MOVES}; got {move!r}.")
        if rule in _DISCRETE_ONLY_RULES:
            raise NotImplementedError(
                f"rule={rule!r}: the heat-bath update for a continuous state "
                "is a conditional resample from the local von-Mises "
                "distribution, which is not wired. Use rule='glauber' (the "
                "Barker acceptance is valid for any symmetric proposal) or "
                "rule='metropolis'."
            )
        if upd_mode == "gillespie":
            raise NotImplementedError(
                "upd_mode='gillespie': the shared BKL engine (ref M3) "
                "assumes a deterministic local proposal, so it is not "
                "available for a continuous state (an angular rate density "
                "is not wired)."
            )
        if tie_flip_p is not None:
            raise ValueError(
                "tie_flip_p has no meaning for a continuous state: "
                "ΔE = 0 has measure zero under the angular proposal."
            )
        # Axis interactions (mirrors PottsMetropolis, plan §9): cluster moves
        # carry their own intrinsic schedule/acceptance, so the inapplicable
        # single-site axes are rejected — never silently ignored.
        if move in ("wolff", "sw"):
            if upd_mode != XY_UPD_MODE_DEFAULT:
                raise ValueError(
                    f"upd_mode does not apply to move={move!r} (the move "
                    "carries its own schedule); leave upd_mode at "
                    f"{XY_UPD_MODE_DEFAULT!r}."
                )
            if order != XY_ORDER_DEFAULT:
                raise ValueError(
                    f"order does not apply to move={move!r}; leave order at "
                    f"{XY_ORDER_DEFAULT!r}."
                )
            if rule != XY_RULE_DEFAULT:
                raise ValueError(
                    f"rule does not apply to move={move!r}: cluster moves "
                    "define their own (Fortuin-Kasteleyn) acceptance."
                )
        elif move == "exchange" and upd_mode != XY_UPD_MODE_DEFAULT:
            raise ValueError(
                "move='exchange' has asynchronous kinetics only (the "
                "swap-attempt sequence is its schedule); leave upd_mode "
                f"at {XY_UPD_MODE_DEFAULT!r}."
            )
        # The proposal width drives move='single' only; an explicit value on
        # any other move would be silently unused — reject it instead.
        if move != "single" and delta is not None:
            raise ValueError(
                f"delta (the single-site proposal width) does not apply to "
                f"move={move!r}."
            )

        self.rule = rule
        self.upd_mode = upd_mode
        self.order = order
        self.move = move
        self.tie_flip_p = None
        self.cluster_sizes: list[int] | None = None
        self.cluster_counts: list[int] | None = None

        super().__init__(
            sg,
            T=T,
            delta=XY_DELTA_DEFAULT if delta is None else float(delta),
            **kwargs,
        )

        # Setup-time capability check: constructing the engine validates
        # upd_mode/order/rule/T combinations NOW, not mid-run. The cluster
        # engines additionally need the couplings, so for those moves the
        # same validation happens in _make_engine on first use, still before
        # any sweep runs.
        if move == "single":
            self._make_engine()

    def _make_engine(self):
        if self.move == "single":
            return ThermalEngine(
                self,
                rule=self.rule,
                upd_mode=self.upd_mode,
                T=self.T,
                order=self.order,
                tie_flip_p=None,
            )
        self._ensure_couplings()
        if self.move in ("wolff", "sw"):
            # Embedded-Ising cluster moves (M6): one reflection direction
            # per cluster (wolff, redrawn after each flip) or per sweep
            # (sw, drawn in edge_satisfactions). Signed couplings are fine.
            self._redraw_dir_on_flip = self.move == "wolff"
            self._draw_reflection()
            if self.move == "wolff":
                return WolffEngine(self, T=self.T)
            return SwendsenWangEngine(self, T=self.T)
        # exchange
        return KawasakiEngine(
            self,
            rule=self.rule,
            T=self.T,
            order=self.order,
            tie_flip_p=None,
        )

    # ------------------------------------------------------------------
    # Native pybind backend (xy_sampling kernel; guards below)
    # ------------------------------------------------------------------
    def _pb_check_supported(self) -> None:
        """What the compiled kernel (``LRGSG_xy.c``) cannot represent is a
        hard capability error (invariant #3): it implements the Metropolis
        acceptance only and visits sites uniformly at random WITH
        replacement (no other order)."""
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
        if XY_OBS_SNAPSHOTS in self._selected_obs:
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
        mod = load_xy_native()
        ni, nw, nptr = self._pb_csr()
        e0 = self.energy_intensive()
        m0 = self.order_parameter()
        theta_out, ene, magn = mod.xy_sampling(
            np.ascontiguousarray(self.s, dtype=np.float64),
            ni,
            nw,
            nptr,
            float(self.T),
            float(self.delta),
            int(self.steps),
            int(self.seed),
            True,
        )
        self.s = np.asarray(theta_out, dtype=np.float64)
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
