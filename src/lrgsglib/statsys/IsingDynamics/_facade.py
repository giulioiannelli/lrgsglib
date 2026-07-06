"""Strangler facade: legacy ``IsingDynamics(...)`` vocabulary → scheme classes.

``ising_from_legacy(sg, ...)`` accepts the legacy god-object's constructor
vocabulary — ``runlang`` codes (``'C0E'``/``'py'``/``'pb_met'``/
``'py_topo_cem'``/...), ``save`` letters (E/S/K/V/H), ``eqSTEP``/
``nstepsIsing``, ``sa_enabled``/``pt_enabled`` and the ``sa_*``/``pt_*``/
``topo_*``/``cem_*`` parameter families — emits ONE ``DeprecationWarning``
naming the translation, and returns the configured NEW scheme class
(plan §10 endpoint table, §11 migration).

Strangler-pattern boundary (deliberate):

- The legacy class itself stays untouched and fully working — this factory
  is ADDITIVE, the new front door for old call sites. Gutting
  ``IsingDynamics.__init__`` into a shim is a later, destructive step.
- C-subprocess codes (``'C...'``) are NOT translated: the C backend stays
  legacy-only (decision D-B5, frozen argv) — a clear error points back to
  the legacy class.
- Anything the new classes cannot represent (save letters K/V/H,
  ``topo_proposal_mode != 'weighted'``) is a hard capability error, never a
  silent drop (invariant #3).

Legacy TFCA is reproduced exactly: a constant-``T`` quench schedule
(``ISING_TFCA_QUENCH_T`` for ``n_temperatures * steps_per_T`` sweeps) under
the spectral field.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Any

import numpy as np

from .defaults import (
    ISING_OBS_ENERGY,
    ISING_OBS_MAGN,
    ISING_OBS_SNAPSHOTS,
    ISING_TFCA_QUENCH_T,
)
from .IsingBase import IsingBase
from .IsingCEM import IsingCEM
from .IsingDynamics import _resolve_runlang
from .IsingMetropolis import IsingMetropolis
from .IsingParallelTempering import IsingParallelTempering
from .IsingSimulatedAnnealing import IsingSimulatedAnnealing

if TYPE_CHECKING:
    from ...graphs.protocols import DynamicsGraphProtocol as SignedGraph

__all__ = ["ising_from_legacy"]

#: Legacy save letter -> new observable names. K (cluster magnetization),
#: V (eigenvector overlap) and H (field file) have no new-class observable
#: yet and raise; cluster statistics live on the cluster moves' attributes.
_SAVE_LETTER_OBSERVABLES: dict[str, tuple[str, ...]] = {
    "E": (ISING_OBS_ENERGY, ISING_OBS_MAGN),
    "S": (ISING_OBS_SNAPSHOTS,),
}

#: Legacy ``upd_mode`` strings -> the async ``order`` axis.
_UPD_MODE_ORDER: dict[str, str] = {
    "asynchronous": "random",
    "sequential": "typewriter",
}


def _resolve_save_letters(save: str | None) -> tuple[str, ...] | None:
    """Map the legacy ``save`` letter string to ``observables=`` names."""
    if save is None:
        return None
    observables: list[str] = []
    for letter in save.strip().upper():
        try:
            names = _SAVE_LETTER_OBSERVABLES[letter]
        except KeyError:
            raise NotImplementedError(
                f"save letter {letter!r} has no new-class observable: "
                "K (cluster magnetization) is recorded on the cluster "
                "moves' attributes, V (eigenvector overlap) and H (field "
                "file) are not wired; supported letters: "
                f"{tuple(_SAVE_LETTER_OBSERVABLES)}."
            ) from None
        for name in names:
            if name not in observables:
                observables.append(name)
    return tuple(observables)


def ising_from_legacy(
    sg: "SignedGraph",
    T: float = 0.0,
    *,
    runlang: str = "py",
    save: str | None = None,
    steps: int | None = None,
    nstepsIsing: int | None = None,
    eqSTEP: int | None = None,
    upd_mode: str = "asynchronous",
    # Simulated Annealing family
    sa_enabled: bool = False,
    T_init: float | None = None,
    T_final: float | None = None,
    cooling_schedule: str | None = None,
    cooling_rate: float | None = None,
    steps_per_T: int | None = None,
    n_temperatures: int | None = None,
    custom_T_schedule: np.ndarray | None = None,
    # Parallel Tempering family
    pt_enabled: bool = False,
    n_replicas: int | None = None,
    T_min: float | None = None,
    T_max: float | None = None,
    T_ladder_type: str | None = None,
    custom_T_ladder: np.ndarray | None = None,
    steps_per_exchange: int | None = None,
    n_exchanges: int | None = None,
    # Topological (spectral) family
    topo_n_modes: int | None = None,
    topo_sigma_init: float | None = None,
    topo_proposal_mode: str | None = None,
    topo_chunk_size: int | None = None,
    topo_polish: bool | None = None,
    topo_polish_sweeps: int | None = None,
    topo_tau: float | None = None,
    topo_field_strength: float | None = None,
    # Cross-Entropy Method family
    cem_iter: int | None = None,
    cem_pop_size: int | None = None,
    cem_elite_frac: float | None = None,
    cem_init_sigma: float | None = None,
    cem_smoothing: float | None = None,
    cem_sigma_floor: float | None = None,
    cem_sigma_ceiling: float | None = None,
    cem_restarts: int | None = None,
    cem_greedy: bool | None = None,
    cem_greedy_sweeps: int | None = None,
    **kwargs: Any,
) -> IsingBase:
    """Translate a legacy ``IsingDynamics(...)`` call into a new scheme class.

    Every parameter keeps its legacy name and meaning; ``None`` means "not
    passed" and leaves the new class's default in force (the defaults were
    lifted from the legacy constructor, so the two agree). Extra ``kwargs``
    (``ic``, ``seed``, ``savedisk``, ``coupling_norm``, ...) are forwarded.

    Returns the CONFIGURED instance — call ``run()`` on it as usual; the
    legacy ``init_ising_dynamics()`` pre-step is no longer needed.
    """
    rl = _resolve_runlang(runlang)
    if rl.startswith("c_"):
        raise NotImplementedError(
            f"runlang={runlang!r} is a C-subprocess code: the C backend "
            "stays on the legacy IsingDynamics class (frozen interface, "
            "decision D-B5) — keep constructing IsingDynamics directly for "
            "C runs."
        )
    # 'py_topo_met' -> ('py', 'topo_met'); an unprefixed algorithm token
    # ('topo_met', 'wolff') ran on the python backend in the legacy
    # dispatcher, so python is the default here too.
    head, _, tail = rl.partition("_")
    if head in ("py", "pb", "cu"):
        backend, algo = head, (tail or "met")
    else:
        backend, algo = "py", rl

    try:
        order = _UPD_MODE_ORDER[upd_mode]
    except KeyError:
        raise ValueError(
            f"Unknown legacy upd_mode {upd_mode!r}; expected one of "
            f"{tuple(_UPD_MODE_ORDER)}."
        ) from None

    resolved_steps = steps if steps is not None else nstepsIsing
    resolved_steps = resolved_steps if resolved_steps is not None else eqSTEP

    observables = _resolve_save_letters(save)

    common: dict[str, Any] = dict(kwargs)
    common["runlang"] = backend
    if observables is not None:
        common["observables"] = observables

    def _put(target: dict, **maybe: Any) -> dict:
        """Copy the non-None entries (unset legacy params keep the new
        class's defaults, which were lifted from the legacy ctor)."""
        target.update({k: v for k, v in maybe.items() if v is not None})
        return target

    # The scheme is named by the runlang's algorithm token; the legacy
    # sa_enabled/pt_enabled flags only upgrade a plain Metropolis code
    # (mirrors the old run() precedence: an explicit algorithm always wins).
    if algo == "met":
        if pt_enabled:
            algo = "pt"
        elif sa_enabled:
            algo = "sa"

    match algo:
        case "topo_cem":
            target_cls: type[IsingBase] = IsingCEM
            new_kwargs = _put(
                common,
                n_modes=topo_n_modes,
                n_iter=cem_iter,
                pop_size=cem_pop_size,
                elite_frac=cem_elite_frac,
                init_sigma=cem_init_sigma,
                smoothing=cem_smoothing,
                sigma_floor=cem_sigma_floor,
                sigma_ceiling=cem_sigma_ceiling,
                restarts=cem_restarts,
                greedy=cem_greedy,
                greedy_sweeps=cem_greedy_sweeps,
            )
        case "topo_met":
            if topo_proposal_mode not in (None, "weighted"):
                raise NotImplementedError(
                    f"topo_proposal_mode={topo_proposal_mode!r}: the "
                    "spectral move proposes RBIM-energy-weighted modes only."
                )
            target_cls = IsingMetropolis
            new_kwargs = _put(
                common,
                T=T,
                move="spectral",
                steps=resolved_steps,
                spectral_n_modes=topo_n_modes,
                spectral_sigma_init=topo_sigma_init,
                spectral_chunk_size=topo_chunk_size,
                spectral_polish=topo_polish,
                spectral_polish_sweeps=topo_polish_sweeps,
            )
        case "topo_fca":
            # Legacy TFCA: sudden quench at a constant near-zero
            # temperature under the spectral field, for
            # n_temperatures * steps_per_T sweeps.
            n_stages = n_temperatures if n_temperatures is not None else 100
            target_cls = IsingSimulatedAnnealing
            new_kwargs = _put(
                common,
                T_schedule=np.full(int(n_stages), ISING_TFCA_QUENCH_T),
                steps_per_T=steps_per_T,
                field="spectral",
                order=order,
                spectral_n_modes=topo_n_modes,
                spectral_tau=topo_tau,
                spectral_field_strength=topo_field_strength,
            )
        case "wolff" | "sw":
            target_cls = IsingMetropolis
            new_kwargs = _put(common, T=T, move=algo, steps=resolved_steps)
        case "pt":
            target_cls = IsingParallelTempering
            new_kwargs = _put(
                common,
                n_replicas=n_replicas,
                T_min=T_min,
                T_max=T_max,
                ladder=T_ladder_type,
                T_ladder=custom_T_ladder,
                steps_per_exchange=steps_per_exchange,
                n_exchanges=n_exchanges,
                order=order,
            )
        case "sa":
            target_cls = IsingSimulatedAnnealing
            new_kwargs = _put(
                common,
                T_init=T_init,
                T_final=T_final,
                schedule=cooling_schedule,
                cooling_rate=cooling_rate,
                steps_per_T=steps_per_T,
                n_temperatures=n_temperatures,
                T_schedule=custom_T_schedule,
                order=order,
            )
        case "met":
            target_cls = IsingMetropolis
            new_kwargs = _put(common, T=T, steps=resolved_steps, order=order)
        case _:
            # No silent Metropolis fallback for a typo'd algorithm token.
            raise ValueError(
                f"Unrecognized legacy algorithm {algo!r} in runlang="
                f"{runlang!r} (resolved {rl!r})."
            )

    warnings.warn(
        f"The legacy IsingDynamics vocabulary is deprecated; this call maps "
        f"to {target_cls.__name__}(...). Construct it directly (plan §10 "
        "endpoint table).",
        DeprecationWarning,
        stacklevel=2,
    )
    return target_cls(sg, **new_kwargs)
