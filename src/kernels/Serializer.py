"""Shared helpers for serializer-style launcher scripts."""

from __future__ import annotations

from typing import Callable, Iterable, Sequence

import numpy as np

__all__ = ["build_memory_function", "build_jobname"]


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
