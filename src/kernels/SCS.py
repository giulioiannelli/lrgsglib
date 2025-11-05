"""Helper utilities for simulations on generalized SCS networks."""

from __future__ import annotations

import ast
from dataclasses import dataclass
from typing import Any

import numpy as np

from lrgsglib import SCSGeneralizedNN

__all__ = [
    "SCSParameters",
    "build_scs_kwargs",
    "normalize_diagonal_argument",
    "resolve_backend",
    "resolve_float_type",
    "probe_output_graph",
]


@dataclass(frozen=True)
class SCSParameters:
    """Container gathering the user supplied SCS model parameters."""

    N: int
    gamma: float
    J: float
    g: float
    diagonal: Any
    workdir: str | None
    seed: int | None


_VALID_BACKENDS = {"numpy", "cupy"}
_VALID_FLOAT_TYPES = {"float32": np.float32, "float64": np.float64}


def resolve_backend(backend: str) -> str:
    """Validate and normalise the backend identifier."""

    backend_norm = backend.lower()
    if backend_norm not in _VALID_BACKENDS:
        raise ValueError(f"Unsupported backend '{backend}'. Expected one of {_VALID_BACKENDS}.")
    if backend_norm == "cupy":
        try:
            import cupy  # noqa: F401
        except ImportError as exc:  # pragma: no cover - informative failure path
            raise RuntimeError(
                "Backend 'cupy' requested but CuPy is not available. Install cupy or choose 'numpy'."
            ) from exc
    return backend_norm


def resolve_float_type(float_type: str) -> type:
    """Map a float type string to the corresponding NumPy dtype."""

    try:
        return _VALID_FLOAT_TYPES[float_type]
    except KeyError as exc:
        raise ValueError(
            f"Unsupported float type '{float_type}'. Expected one of {tuple(_VALID_FLOAT_TYPES)}."
        ) from exc


def normalize_diagonal_argument(value: Any, *, size: int | None = None) -> float | np.ndarray | None:
    """Normalise the diagonal argument for :class:`SCSGeneralizedNN` construction."""

    if value is None:
        return None

    if isinstance(value, (int, float, np.floating)):
        return float(value)

    if isinstance(value, np.ndarray):
        if value.ndim != 1:
            raise ValueError("Diagonal array must be one-dimensional.")
        if size is not None and value.size != size:
            raise ValueError(f"Diagonal array must contain exactly {size} entries.")
        return value.astype(float, copy=False)

    if isinstance(value, (list, tuple)):
        arr = np.asarray(value, dtype=float)
        if size is not None and arr.size != size:
            raise ValueError(f"Diagonal array must contain exactly {size} entries.")
        return arr

    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return None
        lowered = stripped.lower()
        if lowered in {"none", "null", "nan"}:
            return None
        try:
            literal = ast.literal_eval(stripped)
        except (ValueError, SyntaxError):
            literal = stripped
        return normalize_diagonal_argument(literal, size=size)

    raise TypeError("Diagonal specification must be a scalar, iterable, or string literal.")


def build_scs_kwargs(
    params: SCSParameters,
    *,
    J0_value: float,
    make_dir_tree: bool = True,
    seed: int | None = None,
) -> dict[str, Any]:
    """Compose keyword arguments for :class:`SCSGeneralizedNN` construction."""

    diagonal = normalize_diagonal_argument(params.diagonal, size=params.N)
    kwargs: dict[str, Any] = {
        "N": params.N,
        "gamma": params.gamma,
        "J0": float(J0_value),
        "J": params.J,
        "g": params.g,
        "diagonal": diagonal,
        "make_dir_tree": make_dir_tree,
    }
    if params.workdir:
        kwargs["sgpathn"] = params.workdir
    effective_seed = params.seed if seed is None else seed
    if effective_seed is not None:
        kwargs["seed"] = effective_seed
    return kwargs


def probe_output_graph(
    params: SCSParameters,
    *,
    base_J0: float,
) -> SCSGeneralizedNN:
    """Instantiate a minimal SCS graph to initialise filesystem paths."""

    kwargs = build_scs_kwargs(params, J0_value=base_J0, make_dir_tree=True)
    return SCSGeneralizedNN(only_const_mode=True, **kwargs)
