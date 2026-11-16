"""Exact uniform sampler of simple k-regular graphs (pairing model).

Each node gets ``k`` stubs; a uniform random perfect matching of the
``n*k`` stubs is drawn, and the whole draw is rejected if it contains a
self-loop or a multi-edge. Every simple k-regular graph arises from
exactly ``(k!)**n`` matchings, so an accepted draw is uniform over the
labelled simple k-regular graphs on ``n`` nodes.

The acceptance rate falls fast with ``k`` (``expected_pairing_draws``), so
this sampler is meant for small ``k``. Pure numpy: the same seed gives the
same edge array on every engine. The graph classes draw it from
``pairing_stream(seed)``, a child of the seed, so the topology stays
independent of the flip RNG they seed with the bare ``seed``.
"""

from __future__ import annotations

import math
import sys
from typing import Optional, Union

import numpy as np
from numpy.typing import NDArray

from ...config.const import KREG_PAIRING_MAX_ATTEMPTS, KREG_PAIRING_STREAM

__all__ = [
    "canonical_simple_edges",
    "expected_pairing_draws",
    "pairing_stream",
    "random_regular_edges",
    "validate_regular_params",
]

SeedLike = Union[None, int, np.random.SeedSequence, np.random.Generator]


def validate_regular_params(n: int, k: int) -> None:
    """Raise ValueError unless a simple k-regular graph on n nodes exists."""
    if k < 0 or n < 1:
        raise ValueError(f"need n >= 1 and k >= 0, got n={n}, k={k}")
    if k >= n:
        raise ValueError(f"k must be less than n, got k={k}, n={n}")
    if (n * k) % 2 != 0:
        raise ValueError(f"n*k must be even, got n={n}, k={k}, n*k={n*k}")


def pairing_stream(seed: Optional[int]) -> Optional[np.random.SeedSequence]:
    """Seed of the pairing draw of a graph built with ``seed`` (None stays).

    The child ``SeedSequence(seed, spawn_key=(KREG_PAIRING_STREAM,))`` is
    independent of ``default_rng(seed)``, the stream the graph classes use
    for their edge flips.
    """
    if seed is None:
        return None
    return np.random.SeedSequence(seed, spawn_key=(KREG_PAIRING_STREAM,))


_LOG_FLOAT_MAX = math.log(sys.float_info.max)


def expected_pairing_draws(k: int) -> float:
    """Large-n mean number of pairing draws per simple k-regular graph.

    A draw is simple with probability ``exp(-(k**2 - 1) / 4)`` as n grows
    (the classical configuration-model limit; the measured means in
    ``config/const.py`` agree with it), so the mean is the inverse. It is
    ``inf`` once the exponential overflows (k >= 54).
    """
    exponent = (k * k - 1) / 4
    if exponent >= _LOG_FLOAT_MAX:
        return math.inf
    return math.exp(exponent)


def random_regular_edges(
    n: int,
    k: int,
    seed: SeedLike = None,
    max_attempts: int = KREG_PAIRING_MAX_ATTEMPTS,
) -> NDArray[np.int64]:
    """Draw a uniform simple k-regular graph on ``n`` nodes.

    Parameters
    ----------
    n, k : int
        Number of nodes and common degree (``n*k`` even, ``k < n``).
    seed : int, SeedSequence, Generator or None
        Anything ``numpy.random.default_rng`` accepts. A Generator is
        consumed in place.
    max_attempts : int
        Whole draws tried before giving up.

    Returns
    -------
    ndarray of int64, shape (n*k/2, 2)
        Edges with ``u < v`` in each row, rows sorted lexicographically.

    Raises
    ------
    ValueError
        If no simple k-regular graph on ``n`` nodes exists.
    RuntimeError
        If no simple draw was found within ``max_attempts``.
    """
    validate_regular_params(n, k)
    rng = np.random.default_rng(seed)
    stubs = np.repeat(np.arange(n, dtype=np.int64), k)
    for _ in range(max_attempts):
        edges = _simple_matching(stubs, n, rng)
        if edges is not None:
            return edges
    raise RuntimeError(
        f"no simple {k}-regular graph on {n} nodes in {max_attempts} draws"
    )


def canonical_simple_edges(
    pairs: NDArray[np.integer], n: int
) -> Optional[NDArray[np.int64]]:
    """Orient ``(m, 2)`` node pairs as ``u < v`` and sort the rows.

    Returns None if ``pairs`` has a self-loop or a repeated edge. Raises
    ``ValueError`` on endpoints outside ``0..n-1``.
    """
    pairs = np.asarray(pairs, dtype=np.int64).reshape(-1, 2)
    if pairs.size and (pairs.min() < 0 or pairs.max() >= n):
        raise ValueError(f"edge endpoints must lie in 0..{n - 1}")
    u, v = pairs.min(axis=1), pairs.max(axis=1)
    if np.any(u == v):
        return None
    key = u * n + v
    order = np.argsort(key, kind="stable")
    key = key[order]
    if np.any(key[1:] == key[:-1]):
        return None
    return np.column_stack((u[order], v[order]))


def _simple_matching(
    stubs: NDArray[np.int64], n: int, rng: np.random.Generator
) -> Optional[NDArray[np.int64]]:
    """One uniform stub matching as sorted edges, or None if not simple."""
    return canonical_simple_edges(rng.permutation(stubs).reshape(-1, 2), n)
