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
