"""Solver-backend registry for dynamics models.

The dynamics-side analogue of ``graphs/_engine.py``. Maps
``(model_name, SolverBackend) -> Solver`` (or a lazy factory returning one), so a
model's ``run()`` resolves its backend through one shared mechanism instead of
re-implementing a ``runlang`` if/elif chain. ``list_solvers_for_model`` gives the
capability table a UI can render; ``is_backend_available`` probes whether a
backend's runtime (compiled ``.so`` / C binary / CuPy) is actually present.

Example
-------
>>> from lrgsglib.statsys._solver import SolverBackend
>>> from lrgsglib.statsys._solver_engine import get_solver, list_solvers_for_model
>>> list_solvers_for_model("voter")
[<SolverBackend.PY: 'py'>, <SolverBackend.PB: 'pb'>, ...]
>>> solver = get_solver("voter", SolverBackend.PB)
"""

from __future__ import annotations

import os
import warnings
from typing import Callable, Optional, Union

from ._solver import Solver, SolverBackend, coerce_backend

__all__ = [
    "register_solver",
    "get_solver",
    "list_solvers_for_model",
    "list_solver_models",
    "get_backend_from_env",
    "is_backend_available",
]

# model_name -> {backend -> Solver instance | zero-arg factory returning one}
_solvers: dict[str, dict[SolverBackend, Union[Solver, Callable[[], Solver]]]] = {}
# Cache of resolved (instantiated) solvers.
_resolved: dict[str, dict[SolverBackend, Solver]] = {}


def register_solver(
    model: str,
    backend: Union[str, SolverBackend],
    solver_or_factory: Union[Solver, Callable[[], Solver]],
) -> None:
    """Register a solver for ``(model, backend)``.

    ``solver_or_factory`` may be a ready ``Solver`` instance or any callable
    (class or factory) returning one — the latter defers construction until the
    backend is first requested (lazy import of a native module, say).
    """
    model = model.lower()
    backend = coerce_backend(backend)
    _solvers.setdefault(model, {})[backend] = solver_or_factory


def _resolve_solver(model: str, backend: SolverBackend) -> Solver:
    cached = _resolved.get(model, {}).get(backend)
    if cached is not None:
        return cached
    obj = _solvers[model][backend]
    if callable(obj):  # a class or factory — instantiate; instances are not callable
        obj = obj()
    _resolved.setdefault(model, {})[backend] = obj
    return obj


def get_solver(
    model: str,
    backend: Optional[Union[str, SolverBackend]] = None,
    fallback: bool = False,
) -> Solver:
    """Resolve the ``Solver`` for ``(model, backend)``.

    Parameters
    ----------
    model : str
        Registered model name (case-insensitive), e.g. ``'voter'``.
    backend : str, SolverBackend, or None
        Desired backend. If ``None``, uses ``LRGSG_SOLVER_BACKEND`` if set, else
        ``SolverBackend.PY``.
    fallback : bool, default False
        If True and the requested backend is not registered for the model, fall
        back to the Python backend (with a warning). If False (default),
        an unregistered ``(model, backend)`` raises rather than silently
        mislabeling the run.
    """
    model = model.lower()
    if backend is None:
        backend = get_backend_from_env() or SolverBackend.PY
    backend = coerce_backend(backend)

    if model not in _solvers:
        raise ValueError(
            f"No solvers registered for model '{model}'. "
            f"Registered models: {list_solver_models()}."
        )
    impls = _solvers[model]
    if backend in impls:
        return _resolve_solver(model, backend)

    if fallback and SolverBackend.PY in impls:
        warnings.warn(
            f"Model '{model}' has no '{backend.value}' solver; "
            f"falling back to 'py'.",
            stacklevel=2,
        )
        return _resolve_solver(model, SolverBackend.PY)

    available = [b.value for b in impls]
    raise ValueError(
        f"Model '{model}' has no '{backend.value}' solver. Available: {available}."
    )


def list_solvers_for_model(model: str) -> list[SolverBackend]:
    """Backends registered for ``model`` (capability table; ignores runtime)."""
    return list(_solvers.get(model.lower(), {}).keys())


def list_solver_models() -> list[str]:
    """All model names with at least one registered solver."""
    return list(_solvers.keys())


def get_backend_from_env() -> Optional[SolverBackend]:
    """Default backend from ``LRGSG_SOLVER_BACKEND`` if set and valid, else None."""
    raw = os.environ.get("LRGSG_SOLVER_BACKEND")
    if raw:
        try:
            return SolverBackend(raw.lower())
        except ValueError:
            warnings.warn(
                f"Invalid LRGSG_SOLVER_BACKEND value: {raw}. "
                f"Valid options: {[b.value for b in SolverBackend]}. Ignoring.",
                stacklevel=2,
            )
    return None


def is_backend_available(model: str, backend: Union[str, SolverBackend]) -> bool:
    """Whether ``(model, backend)`` is registered *and* its runtime is present.

    Delegates to the solver's optional ``is_available()`` (compiled module / C
    binary / CuPy probe); a solver without one is assumed always-available.
    """
    model = model.lower()
    backend = coerce_backend(backend)
    if backend not in _solvers.get(model, {}):
        return False
    solver = _resolve_solver(model, backend)
    check = getattr(solver, "is_available", None)
    if callable(check):
        try:
            return bool(check())
        except Exception:
            return False
    return True
