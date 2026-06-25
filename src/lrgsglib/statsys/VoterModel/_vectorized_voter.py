"""Vectorized synchronous voter sweeps (NumPy / CuPy).

The **synchronous** voter is embarrassingly parallel: every node updates from the
*same* frozen snapshot, so a whole sweep is a handful of CSR gathers / segment
reductions with no inter-node dependence (unlike the asynchronous schedule, and
unlike Ising which needs graph-colouring). One array routine therefore serves
both backends -- pass ``xp=numpy`` (``np``) or ``xp=cupy`` (``cu``) -- across the
full rule family:

* ``linear``    -- each node copies one uniformly chosen signed neighbour.
* ``majority``  -- each node adopts ``sign(Σ_j sign(w_ij) s_j)`` (random tie-break).
* ``qvoter``    -- sample ``q`` neighbours with replacement; copy if signed-unanimous,
  else flip with probability ``eps``.
* ``nonlinear`` -- adopt ``+1`` with probability ``f₊^α/(f₊^α+f₋^α)``, ``f₊`` the
  signed-neighbour ``+1`` fraction (Yang et al. 2012; see ``defaults.py``).

Only the synchronous schedule is vectorized; the asynchronous / link / gillespie
schedules are intrinsically sequential and stay on the py / C / pybind backends.

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
    rule: str = "linear",
    q: int = 2,
    eps: float = 0.0,
    alpha: float = 1.0,
):
    """Run the synchronous voter on CSR adjacency with array module ``xp``.

    Parameters mirror the CSR layout from ``statsys/_csr.py``
    (``indices``/``weights`` length ``2|E|``, ``ptr`` length ``N+1``). ``rule``
    selects the local update (``linear`` / ``majority`` / ``qvoter`` /
    ``nonlinear``); ``q``/``eps`` parametrise the q-voter and ``alpha`` the
    nonlinear rule (ignored by the others).

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
    # source node of each half-edge: needed for the majority / nonlinear segment
    # reductions and for vectorized frustration counting. Built on the host from
    # `ptr` -- cupy.repeat rejects a device array as `repeats`.
    deg_host = np.diff(np.asarray(ptr, dtype=np.int64))
    src = xp.asarray(np.repeat(np.arange(N, dtype=np.int64), deg_host))

    rng = xp.random.default_rng(seed)
    magn_dev = xp.empty(n_sweeps, dtype=xp.float64) if savemagn else None
    nrec = 0
    absorbed_at = None

    def _sweep(state):
        """One synchronous sweep under ``rule``; returns the new state (int8)."""
        if rule == "linear":
            off = (rng.random(N) * degf).astype(xp.int64)     # uniform in [0, deg)
            pick = xp.where(has, starts + off, 0)
            s_new = (sgn[pick] * state[idx[pick]]).astype(xp.int8)
            return xp.where(has, s_new, state)
        # signed neighbour opinion per half-edge (±1), shared by majority/nonlinear
        if rule in ("majority", "nonlinear"):
            op = sgn.astype(xp.int64) * state[idx].astype(xp.int64)
            if rule == "majority":
                h = xp.bincount(src, weights=op.astype(xp.float64), minlength=N)
                new = xp.sign(h)                              # -1 / 0 / +1
                tie = new == 0
                rnd = xp.where(rng.random(N) < 0.5, 1.0, -1.0)
                new = xp.where(tie, rnd, new).astype(xp.int8)
                return xp.where(has, new, state)
            # nonlinear: P(+1) = f₊^α / (f₊^α + f₋^α), f₊ = +1 neighbour fraction
            nplus = xp.bincount(src, weights=(op > 0).astype(xp.float64),
                                minlength=N)
            fp = xp.where(has, nplus / xp.where(has, degf, 1.0), 0.5)
            fm = 1.0 - fp
            num = fp ** alpha
            den = num + fm ** alpha
            pplus = xp.where(den > 0.0, num / xp.where(den > 0.0, den, 1.0), 0.5)
            new = xp.where(rng.random(N) < pplus, 1, -1).astype(xp.int8)
            return xp.where(has, new, state)
        if rule == "qvoter":
            # sample q neighbours WITH replacement; copy if signed-unanimous, else
            # flip with probability eps (q=1 reduces to the linear voter).
            ops = xp.empty((q, N), dtype=xp.int8)
            for k in range(q):
                off = (rng.random(N) * degf).astype(xp.int64)
                pick = xp.where(has, starts + off, 0)
                ops[k] = (sgn[pick] * state[idx[pick]]).astype(xp.int8)
            all_plus = (ops == 1).all(axis=0)
            all_minus = (ops == -1).all(axis=0)
            unanimous = all_plus | all_minus
            unval = xp.where(all_plus, 1, -1).astype(xp.int8)
            flip = rng.random(N) < eps
            new = xp.where(unanimous, unval,
                           xp.where(flip, (-state).astype(xp.int8), state))
            return xp.where(has, new, state)
        raise ValueError(f"unsupported vectorized rule={rule!r}")

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
        s = _sweep(s)

    magn: list[float] = []
    if savemagn:
        md = magn_dev[:nrec]
        magn = (_cp.asnumpy(md) if xp is not np else md).tolist()

    final = s if xp is np else _cp.asnumpy(s)
    return np.asarray(final, dtype=np.int8), magn, absorbed_at
