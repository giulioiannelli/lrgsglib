"""Canonical run-directory naming for statsys outputs (Phase C).

Every run writes ONE directory under the model's dynamics path::

    data/<graph>/<model>/N=.../<run dirname>/{ene.bin, m.bin, sout.bin, cfg.json}

and the run dirname is built here from the model's resolved parameters.

Rules (locked with the user, 2026-07-07):

- ``key=value`` for EVERY token — numbers and categoricals alike
  (``rule=glau_move=wolff_upd=gill``), floats ``%.3g``, integer
  sequences joined with ``'-'`` (heterogeneous q: ``q=2-5``).
- Token order is the model's canonical schema: numeric parameters first
  in order of fundamentality (disorder ``p=``, then model physics, then
  run length), categorical axes after, ``lang=`` and ``s=`` always last.
  The schema owner passes tokens already ordered; this module never
  reorders.
- A token whose value equals its declared default is ELIDED, so a plain
  run reads ``p=0.1_T=2.27_ns=5000_lang=pb_s=91731`` and anything
  non-default becomes visible. Tokens declared with ``ALWAYS`` (seed,
  backend, the fundamental physics numerics) are never elided.
- Categorical values are shortened through ``VOCAB`` (``glauber`` ->
  ``glau``); unknown strings pass through unchanged. Inapplicable axes
  are simply not passed (a ``move=wolff`` run has no ``rule=`` token —
  cluster acceptance is baked into the bond probability).

The dirname is for humans and globs; the ``cfg.json`` sidecar in the
run directory is the reproducibility ground truth.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Iterable, Sequence

__all__ = [
    "ALWAYS",
    "VOCAB",
    "format_value",
    "run_dirname",
    "parse_run_dirname",
    "write_cfg_sidecar",
]

#: Sentinel default: a token declared with ``ALWAYS`` is never elided.
ALWAYS: Any = object()

#: Canonical short forms for categorical token values. Values not listed
#: here pass through unchanged (e.g. ``wolff``, ``sw``, ``sync``).
VOCAB: dict[str, str] = {
    # schemes
    "metropolis": "met",
    "simulated_annealing": "sa",
    "parallel_tempering": "pt",
    "cross_entropy": "cem",
    # acceptance rules
    "glauber": "glau",
    "heatbath": "hb",
    # moves
    "single": "sing",
    "spectral": "spec",
    "exchange": "exch",
    # update modes
    "asynchronous": "async",
    "synchronous": "sync",
    "gillespie": "gill",
    # sweep orders
    "random": "rand",
    "permutation": "perm",
    "typewriter": "typw",
    # initial conditions
    "uniform": "unif",
    # voter rules
    "linear": "lin",
    "majority": "maj",
    "qvoter": "qv",
    "nonlinear": "nlin",
    # schedules / ladders
    "exponential": "exp",
    "geometric": "geo",
    # CP-EI update strategies (the naive default is elided anyway)
    "frontier": "front",
    "adaptive": "adpt",
    # reaction types (ReactionDiffusion)
    "fisher_kpp": "fkpp",
    "bistable": "bist",
    # ODE local terms
    "logistic": "logi",
    "fitzhugh": "fitz",
    # coupling functions
    "product": "prod",
    # frequency / disorder distributions
    "gaussian": "gauss",
    "normal": "gauss",
}


def format_value(value: Any) -> str:
    """Render one token value: floats ``%.3g``, sequences ``'-'``-joined,
    categorical strings shortened through ``VOCAB``."""
    if isinstance(value, bool):
        return str(int(value))
    if isinstance(value, float):
        return f"{value:.3g}"
    if isinstance(value, (tuple, list)):
        return "-".join(format_value(v) for v in value)
    if isinstance(value, str):
        return VOCAB.get(value, value)
    return str(value)


def _values_equal(value: Any, default: Any) -> bool:
    if isinstance(value, (tuple, list)) or isinstance(default, (tuple, list)):
        return tuple(value) == tuple(default)
    return bool(value == default)


def run_dirname(tokens: Iterable[Sequence[Any]]) -> str:
    """Join ``(key, value[, default])`` triples into the run dirname.

    ``value is None`` skips the token (inapplicable axis); a value equal
    to its ``default`` is elided; no third element (or ``ALWAYS``) means
    the token is always written.
    """
    parts: list[str] = []
    for token in tokens:
        key, value = token[0], token[1]
        default = token[2] if len(token) > 2 else ALWAYS
        if value is None:
            continue
        if default is not ALWAYS and _values_equal(value, default):
            continue
        parts.append(f"{key}={format_value(value)}")
    return "_".join(parts)


def parse_run_dirname(name: str) -> dict[str, str]:
    """Inverse of ``run_dirname`` at the string level: ``key -> raw value``.

    Values are returned as the formatted strings (``'2.27'``, ``'2-5'``,
    ``'glau'``); the ``cfg.json`` sidecar holds the typed originals.
    """
    out: dict[str, str] = {}
    for part in name.split("_"):
        if "=" not in part:
            raise ValueError(
                f"malformed run-dirname token {part!r} in {name!r} "
                f"(every token is key=value)"
            )
        key, value = part.split("=", 1)
        out[key] = value
    return out


def _json_default(obj: Any) -> Any:
    if hasattr(obj, "tolist"):  # numpy scalars and arrays
        return obj.tolist()
    if isinstance(obj, Path):
        return str(obj)
    return str(obj)


def write_cfg_sidecar(
    rundir: Path,
    graph: dict[str, Any],
    model: dict[str, Any],
    run: dict[str, Any],
) -> Path:
    """Write the ``cfg.json`` reproducibility sidecar into ``rundir``.

    Three blocks: ``graph`` (substrate identity incl. the disorder-
    realization seed — ``p`` names the ensemble, the seed pins the
    realization), ``model`` (every constructor argument, elided defaults
    included), ``run`` (steps, observables, backend, dynamics seed,
    library version, timestamp). Atomic write (temp + rename).
    """
    payload = {"graph": graph, "model": model, "run": run}
    path = Path(rundir) / "cfg.json"
    tmp = path.with_suffix(".json.tmp")
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True, default=_json_default)
        fh.write("\n")
    os.replace(tmp, path)
    return path
