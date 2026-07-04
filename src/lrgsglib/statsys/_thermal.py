"""Generic finite-temperature MCMC engine (shared spin-model home).

PHASE-0 SKELETON — shapes only, no kernels yet. The acceptance-rule bodies and
the compiled sweep land with the Ising vertical slice (Phase 1), AFTER the MCMC
references are added to ``.agents/ref/00-REFERENCES.md`` (citation debt logged
in the Ising plan §12): every acceptance rule ships with its primary-source
citation or it does not ship.

Why this file exists (plan §4, decision D-B2): the sweep / single-flip /
acceptance machinery is a *spin-model* layer, not a base-class layer. It is
written ONCE here against the substrate contract below and consumed by every
spin model — Ising (``BinDynSys``-based) today, Potts/XY/Heisenberg/
MultiSpecies (``VecDynSys``-based) as they modernize. Flow models (Kuramoto,
reaction-diffusion, coupled ODE) and the walker never import it.

Design contract (plan §5, configure-then-run): everything is resolved ONCE at
setup — the acceptance rule, the update schedule, the move generator — into
plain callables; the hot loop only calls them. No per-step branching on
strings, ever.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any, Protocol, runtime_checkable

import numpy as np

__all__ = [
    "THERMAL_RULES",
    "UPDATE_MODES",
    "ThermalSubstrate",
    "resolve_acceptance",
    "ThermalEngine",
]

#: Acceptance rules the engine will support. Kernel bodies + citations land in
#: Phase 1; keeping the vocabulary here lets scheme classes validate early.
THERMAL_RULES: tuple[str, ...] = ("metropolis", "glauber", "heatbath")

#: Update schedules (plan §9): random-sequential single-site, checkerboard /
#: sublattice-parallel, and rejection-free continuous-time (BKL) kinetics.
UPDATE_MODES: tuple[str, ...] = ("async", "sync", "gillespie")


@runtime_checkable
class ThermalSubstrate(Protocol):
    """The substrate contract (plan §3) a spin model exposes to the engine.

    The engine never sees model internals — only these six operations. This is
    what makes the engine write-once: ``IsingBase`` implements them on an int8
    ±1 state; a future ``PottsBase`` implements them on an integer-label state;
    the engine is unchanged.
    """

    #: Current state array (any dtype/shape the model chooses).
    s: np.ndarray

    def propose_local(self, nd: int) -> Any:
        """Propose a new local configuration at node *nd* (e.g. the flipped
        spin, a new Potts label, a rotated XY angle)."""
        ...

    def delta_E(self, nd: int, proposal: Any) -> float:
        """Energy change if *proposal* were committed at node *nd* (includes
        the external-field term)."""
        ...

    def commit(self, nd: int, proposal: Any) -> None:
        """Apply *proposal* at node *nd* (in place)."""
        ...

    def order_parameter(self) -> float:
        """The model's scalar order parameter for the current state."""
        ...

    def bond_activation(self, u: int, v: int) -> float:
        """Bond-activation probability for cluster moves (Wolff / SW); models
        without cluster moves may raise ``NotImplementedError``."""
        ...


def resolve_acceptance(rule: str, T: float) -> Callable[[float], float]:
    """Resolve an acceptance rule name to a compiled ``ΔE -> p_accept`` callable.

    Called ONCE at setup (never in the hot loop). Phase 1 fills in the rule
    bodies together with their primary-source citations in
    ``.agents/ref/00-REFERENCES.md``.
    """
    if rule not in THERMAL_RULES:
        raise ValueError(
            f"Unknown acceptance rule {rule!r}; expected one of {THERMAL_RULES}."
        )
    raise NotImplementedError(
        "Phase-0 skeleton: acceptance kernels land with the Ising vertical "
        "slice (Phase 1), after the MCMC citation debt is settled."
    )


class ThermalEngine:
    """Configure-then-run MCMC driver over a :class:`ThermalSubstrate`.

    Parameters are resolved to callables in :meth:`compile_step`; the returned
    step function advances one sweep with zero per-step branching. Scheme
    classes (``IsingMetropolis``, ``IsingSimulatedAnnealing``, ...) own the
    temperature *policy* (fixed T, T(t) schedule, replica ladder) and drive
    this engine; the engine owns the *mechanics* of a sweep.
    """

    def __init__(
        self,
        substrate: ThermalSubstrate,
        *,
        rule: str,
        upd_mode: str,
        T: float,
    ) -> None:
        if rule not in THERMAL_RULES:
            raise ValueError(
                f"Unknown acceptance rule {rule!r}; expected one of {THERMAL_RULES}."
            )
        if upd_mode not in UPDATE_MODES:
            raise ValueError(
                f"Unknown upd_mode {upd_mode!r}; expected one of {UPDATE_MODES}."
            )
        self.substrate = substrate
        self.rule = rule
        self.upd_mode = upd_mode
        self.T = float(T)

    def compile_step(self) -> Callable[[], None]:
        """Resolve (rule × upd_mode) ONCE into a zero-argument sweep callable.

        Phase-1 deliverable; declared now so scheme classes can code against
        the engine surface.
        """
        raise NotImplementedError(
            "Phase-0 skeleton: the compiled sweep lands with the Ising "
            "vertical slice (Phase 1)."
        )
