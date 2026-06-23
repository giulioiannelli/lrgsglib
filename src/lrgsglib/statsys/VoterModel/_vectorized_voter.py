"""Vectorized synchronous linear-voter sweeps (NumPy / CuPy).

The **synchronous** linear voter is embarrassingly parallel: every node copies a
uniformly chosen signed neighbour from the *same* frozen snapshot, so a whole
sweep is a single CSR gather with no inter-node dependence (unlike the
asynchronous schedule, and unlike Ising which needs graph-colouring). One array
routine therefore serves both backends -- pass ``xp=numpy`` (``np``) or
``xp=cupy`` (``cu``).

Signed substrate: the effective neighbour opinion is ``sign(w_ij) * s_j`` (a
negative edge is an anti-copy), matching the other VoterModel backends.
"""

from __future__ import annotations

import numpy as np

try:  # optional GPU dependency
    import cupy as _cp  # noqa: F401
    CUPY_AVAILABLE = True
except Exception:  # pragma: no cover - depends on the host
    _cp = None
    CUPY_AVAILABLE = False


def run_vectorized_sync(
    xp,
    s0: np.ndarray,
    indices: np.ndarray,
    weights: np.ndarray,
    ptr: np.ndarray,
    n_sweeps: int,
    seed: int,
    savemagn: bool,
    absorbing_check: bool,
    absorbing_every: int,
):
    """Run the synchronous linear voter on CSR adjacency with array module ``xp``.

    Parameters mirror the CSR layout from ``statsys/_csr.py``
    (``indices``/``weights`` length ``2|E|``, ``ptr`` length ``N+1``).

    Returns ``(final_spins, magn, absorbed_at)`` with ``final_spins`` an
    ``int8`` host ``ndarray``, ``magn`` a Python list (empty unless
    ``savemagn``), and ``absorbed_at`` the sweep index of the first
    zero-frustration configuration or ``None``.
    """
    s = xp.asarray(s0, dtype=xp.int8)
    idx = xp.asarray(indices, dtype=xp.int64)
    w = xp.asarray(weights, dtype=xp.float64)
    sgn = xp.where(w < 0, -1, 1).astype(xp.int8)          # per half-edge sign
    p = xp.asarray(ptr, dtype=xp.int64)

    N = int(s.shape[0])
    deg = p[1:] - p[:-1]
    degf = deg.astype(xp.float64)
    starts = p[:-1]
    has = deg > 0
    # source node of each half-edge (for vectorized frustration counting). Built
    # on the host from `ptr` -- cupy.repeat rejects a device array as `repeats`.
    src = None
    if absorbing_check:
        deg_host = np.diff(np.asarray(ptr, dtype=np.int64))
        src = xp.asarray(np.repeat(np.arange(N, dtype=np.int64), deg_host))

    rng = xp.random.default_rng(seed)
    magn_dev = xp.empty(n_sweeps, dtype=xp.float64) if savemagn else None
    nrec = 0
    absorbed_at = None

    for t in range(n_sweeps):
        if savemagn:
            magn_dev[t] = s.mean()       # stays on device (no per-sweep sync)
            nrec = t + 1
        if absorbing_check and (t % absorbing_every == 0):
            prod = (sgn.astype(xp.int64) * s[src].astype(xp.int64)
                    * s[idx].astype(xp.int64))
            if int((prod < 0).sum()) == 0:   # zero frustrated edges => absorbing
                absorbed_at = t
                break
        # one synchronous sweep: each node copies a random signed neighbour
        off = (rng.random(N) * degf).astype(xp.int64)     # uniform in [0, deg)
        pick = xp.where(has, starts + off, 0)
        s_new = (sgn[pick] * s[idx[pick]]).astype(xp.int8)
        s = xp.where(has, s_new, s)        # isolated nodes unchanged

    magn: list[float] = []
    if savemagn:
        md = magn_dev[:nrec]
        magn = (_cp.asnumpy(md) if xp is not np else md).tolist()

    final = s if xp is np else _cp.asnumpy(s)
    return np.asarray(final, dtype=np.int8), magn, absorbed_at
