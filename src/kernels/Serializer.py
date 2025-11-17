"""Shared helpers for serializer-style launcher scripts."""

from __future__ import annotations

from typing import Callable, Iterable, Sequence

import numpy as np

__all__ = ["build_memory_function", "build_jobname", "_determine_precision", "_format_value_consistently", "_collect_values", "_collect_values_typed"]


def build_memory_function(min_mb: int, max_mb: int, values: Sequence[int]) -> Callable[[int], int]:
    """
    Return an interpolating function that maps parameter values to memory limits.

    When the provided min/max limits are identical, a constant function is returned.
    """
    if min_mb == max_mb or not values:
        return lambda *_: min_mb

    low, high = min(values), max(values)
    if low == high:
        return lambda *_: min_mb

    bounds = [low, high]
    targets = [min_mb, max_mb]

    def _interp(value: int) -> int:
        return int(np.interp(value, bounds, targets))

    return _interp


def build_jobname(
    *,
    program_short: str,
    tokens: Iterable[str],
    job_id: str | None = None,
    delimiter: str = "_",
    prefix_job_id: bool = False,
) -> str:
    """
    Compose a Slurm job name using a program short name and additional tokens.

    Set ``prefix_job_id`` to True if the external job identifier should appear
    before the program short name (matching historical behaviour in some scripts).
    """
    cleaned_tokens = [tok for tok in tokens if tok]
    parts = [program_short, *cleaned_tokens]

    if job_id:
        if prefix_job_id:
            parts = [job_id, *parts]
        else:
            parts = [program_short, job_id, *cleaned_tokens]

    return delimiter.join(parts)


def _determine_precision(values: list[float]) -> int:
    """
    Determine the minimum number of decimal places needed to distinguish all values.
    Returns the number of decimal places needed.
    """
    if len(values) <= 1:
        # For single value or empty, use minimal precision
        return 2

    # Sort values to find minimum differences
    sorted_vals = sorted(set(values))
    if len(sorted_vals) == 1:
        return 2

    # Find the minimum difference between consecutive values
    min_diff = min(abs(sorted_vals[i+1] - sorted_vals[i]) for i in range(len(sorted_vals) - 1))

    # Determine precision needed: we need enough decimals so that min_diff is distinguishable
    if min_diff == 0:
        return 2

    # Calculate how many decimal places we need
    # We want at least 2 significant digits in the minimum difference
    precision = max(2, int(np.ceil(-np.log10(min_diff) + 1)))

    # Cap at reasonable maximum
    return min(precision, 8)


def _format_value_consistently(value: float, precision: int) -> str:
    """
    Format a value with consistent decimal places.
    All values use exactly the same number of decimal places (the precision).
    """
    return f"{value:.{precision}f}"


def _collect_values(
    explicit: Iterable[float] | None,
    linspace_values,
) -> list[float]:
    return _collect_values_typed(explicit, linspace_values, float)


def _collect_values_typed(
    explicit: Iterable[float] | None,
    linspace_values,
    caster: Callable[[float], float],
) -> list[float]:
    values: list[float] = []
    if linspace_values is not None:
        # If linspace is provided, use only linspace values
        values.extend(caster(v) for v in np.asarray(linspace_values, dtype=float))
    elif explicit:
        # Otherwise, use explicit values (which may be defaults from argparse)
        values.extend(caster(v) for v in explicit)
    return values
