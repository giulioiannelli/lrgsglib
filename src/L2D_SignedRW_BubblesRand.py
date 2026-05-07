"""Random-start bubble-touch sweep on a 2D signed lattice.

For each disorder realisation we draw ``n_pairs`` random pairs ``(a_k, b_k)``
uniformly on the torus and run ``n_walkers`` absorbing walkers from each
endpoint.  The two **bubbles** of pair ``k`` are the unions of visited sites
for each family; we record the binary touch flag and the Jaccard overlap

    J(k) = |bubble_A(k) ∩ bubble_B(k)| / |bubble_A(k) ∪ bubble_B(k)|.

Averaging over pairs (and disorder) samples the disordered lattice uniformly:
close pairs almost always touch, antipodal pairs only touch below p_c.  Near
``p_c`` the distribution of `J(k)` across random pair-placements develops the
large fluctuations that mark the percolation-like transition.

Runs through the native ``absorb_walker_sampling`` kernel batched over all
pairs in a single call (`walkers_per_trial = n_walkers` groups the starts into
per-pair trials), so compute cost stays linear in ``n_pairs * n_walkers``
without Python-level loop overhead.  Mirrors the sweep skeleton of
``L2D_SignedRW_Bubbles.py`` (argparse, SIGINT checkpoint, per-(p,d) resume).

Usage
-----
python lrgsglib/src/L2D_SignedRW_BubblesRand.py \
    --side 64 --pflip-list 0.06,0.07,...,0.13 \
    --n-walkers 8192 --n-pairs 410 --n-disorder 4 \
    --workdir bubble_random
"""

from __future__ import annotations

import argparse
import pickle
import signal
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter

import numpy as np


_INTERRUPTED = {"flag": False}


def _handle_stop_signal(signum, frame):  # noqa: ARG001
    _INTERRUPTED["flag"] = True
    print(f"\n[interrupt: signal {signum} received; finishing current p-value "
          "and exiting after checkpoint]", flush=True)


from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.SignedRW._kernel import (
    _kill_masks,
    signed_lattice_tables,
)
from lrgsglib.statsys.SignedRW.SignedWalker import _load_srw_native


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_float_list(s: str) -> list[float]:
    return [float(x) for x in s.split(",") if x]


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--side", "-s", type=int, required=True)
    p.add_argument("--pflip-list", required=True, type=_parse_float_list)
    p.add_argument("--geo", default="sqr")
    p.add_argument("--pbc", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--x-node", choices=["reflect", "absorb"], default="reflect")
    p.add_argument("--n-walkers", "-nw", type=int, default=500,
                   help="Walkers per bubble per random pair.")
    p.add_argument("--n-pairs", "-np", type=int, default=400,
                   help="Random (a,b) pairs sampled per disorder realisation.")
    p.add_argument("--n-disorder", "-nd", type=int, default=4)
    p.add_argument("--coverage", type=float, default=1.0,
                   help="Coverage-stop fraction (1.0 = run until death).")
    p.add_argument("--min-dist", type=int, default=0,
                   help="Optional minimum torus Chebyshev distance between "
                        "start nodes a, b (0 = no constraint).")
    p.add_argument("--n-r-bins", type=int, default=8,
                   help="Number of Chebyshev-distance bins for P_touch(r; p).")
    p.add_argument("--chunk-pairs", type=int, default=0,
                   help="Chunk size for per-call pair batching (0 = all at once). "
                        "Smaller chunks reduce per-call working set and improve "
                        "cache locality for large n_walkers * n_pairs products.")
    p.add_argument("--seed-pairs", type=int, default=2024,
                   help="RNG seed for pair sampling.")
    p.add_argument("--seed-a", type=int, default=42)
    p.add_argument("--seed-b", type=int, default=43)
    p.add_argument("--seed-disorder", type=int, default=42)
    p.add_argument("--time-budget", type=float, default=0.0,
                   help="Wall-clock cap per p value in seconds (0 = disabled).")
    p.add_argument("--workdir", "-wd", default="bubble_random")
    p.add_argument("--tag", default="")
    p.add_argument("--verbose", "-v", action="store_true")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

@dataclass
class _CSRMasks:
    N: int
    neigh_indices: np.ndarray
    neigh_ptr: np.ndarray
    kill_flat: np.ndarray
    refl_flat: np.ndarray


def _csr_masks(lat, x_node: str) -> _CSRMasks:
    nbrs, neg, frust = signed_lattice_tables(lat)
    km, rm = _kill_masks(neg, frust, x_node)
    N, z = nbrs.shape
    return _CSRMasks(
        N=N,
        neigh_indices=np.ascontiguousarray(nbrs.ravel(), dtype=np.int64),
        neigh_ptr=np.arange(0, (N + 1) * z, z, dtype=np.int64),
        kill_flat=np.ascontiguousarray(km.ravel(), dtype=np.uint8),
        refl_flat=np.ascontiguousarray(rm.ravel(), dtype=np.uint8),
    )


def _pair_distances(pairs: np.ndarray, side: int) -> np.ndarray:
    """Torus Chebyshev distance for each (a, b) pair.  Shape (n_pairs,)."""
    ar, ac = pairs[:, 0] // side, pairs[:, 0] % side
    br, bc = pairs[:, 1] // side, pairs[:, 1] % side
    dr = np.minimum(np.abs(ar - br), side - np.abs(ar - br))
    dc = np.minimum(np.abs(ac - bc), side - np.abs(ac - bc))
    return np.maximum(dr, dc)


def _sample_pairs(N: int, side: int, n_pairs: int, min_dist: int,
                  rng: np.random.Generator) -> np.ndarray:
    """Sample `n_pairs` random (a, b) pairs on an L×L torus.

    If ``min_dist > 0``, rejection-sample pairs whose torus Chebyshev
    distance is below the minimum.  Returns an (n_pairs, 2) int64 array.
    """
    if min_dist <= 0:
        a = rng.integers(0, N, size=n_pairs, dtype=np.int64)
        b = rng.integers(0, N, size=n_pairs, dtype=np.int64)
        # Guard against a == b (rare, but degenerate).
        same = a == b
        while same.any():
            b[same] = rng.integers(0, N, size=int(same.sum()), dtype=np.int64)
            same = a == b
        return np.stack([a, b], axis=1)
    pairs = np.empty((n_pairs, 2), dtype=np.int64)
    k = 0
    while k < n_pairs:
        batch = min(4 * n_pairs, 100_000)
        ab = rng.integers(0, N, size=(batch, 2), dtype=np.int64)
        ar, ac = ab[:, 0] // side, ab[:, 0] % side
        br, bc = ab[:, 1] // side, ab[:, 1] % side
        dr = np.minimum(np.abs(ar - br), side - np.abs(ar - br))
        dc = np.minimum(np.abs(ac - bc), side - np.abs(ac - bc))
        dist = np.maximum(dr, dc)
        ok = (dist >= min_dist) & (ab[:, 0] != ab[:, 1])
        idx = np.where(ok)[0]
        take = min(len(idx), n_pairs - k)
        pairs[k:k + take] = ab[idx[:take]]
        k += take
    return pairs


def _bin_by_distance(values: np.ndarray, dists: np.ndarray,
                     bin_edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (mean_per_bin, count_per_bin) for `values` bucketed by `dists`.

    `bin_edges` has length n_bins+1; bins are [edges[i], edges[i+1]).
    Bins with zero count return NaN mean.
    """
    n_bins = len(bin_edges) - 1
    means  = np.full(n_bins, np.nan, dtype=np.float64)
    counts = np.zeros(n_bins, dtype=np.int64)
    idx = np.digitize(dists, bin_edges[1:-1])  # 0..n_bins-1
    for b in range(n_bins):
        sel = idx == b
        if sel.any():
            means[b]  = float(values[sel].mean())
            counts[b] = int(sel.sum())
    return means, counts


def _bubbles_batched(native, csr: _CSRMasks, starts: np.ndarray, *,
                     seed: int, n_walkers: int,
                     stop_thresh: int,
                     chunk_pairs: int = 0) -> np.ndarray:
    """Pybind batched call producing per-pair bubble bitmaps.

    ``starts`` has shape ``(n_pairs * n_walkers,)`` where each contiguous
    block of size ``n_walkers`` shares a start node.  Emits ``(n_pairs, N)``
    uint8 bitmap.

    If ``chunk_pairs > 0`` the pair batch is split into sub-batches of that
    size and the kernel is called once per chunk; the per-call working set
    shrinks from ``n_pairs * n_walkers`` to ``chunk_pairs * n_walkers``, which
    is much more cache-friendly when the product exceeds ~1e7 walker-starts.
    Seeds for sub-batches are offset by the sub-batch index.
    """
    total = starts.size
    if total % n_walkers:
        raise ValueError(
            f"starts.size ({total}) not a multiple of n_walkers ({n_walkers})."
        )
    n_pairs = total // n_walkers
    if chunk_pairs <= 0 or chunk_pairs >= n_pairs:
        out = native.absorb_walker_sampling(
            csr.N,
            csr.neigh_indices, csr.neigh_ptr,
            csr.kill_flat, csr.refl_flat,
            starts,
            int(stop_thresh),
            int(seed),
            int(n_walkers),
            False,
        )
        return out[5]
    bubbles = np.empty((n_pairs, csr.N), dtype=np.uint8)
    base = 0
    sub_idx = 0
    while base < n_pairs:
        end = min(base + chunk_pairs, n_pairs)
        sub_starts = starts[base * n_walkers:end * n_walkers]
        out = native.absorb_walker_sampling(
            csr.N,
            csr.neigh_indices, csr.neigh_ptr,
            csr.kill_flat, csr.refl_flat,
            sub_starts,
            int(stop_thresh),
            int(seed + sub_idx),
            int(n_walkers),
            False,
        )
        bubbles[base:end] = out[5]
        base = end
        sub_idx += 1
    return bubbles


# ---------------------------------------------------------------------------
# Output path
# ---------------------------------------------------------------------------

def _output_path(lat, n_walkers: int, n_pairs: int, tag: str) -> Path:
    base = Path(lat.path_srw) / "bubbles_random"
    base.mkdir(parents=True, exist_ok=True)
    stem = (f"bubble_rand_L={lat.side1}_nw={n_walkers}_np={n_pairs}"
            + (f"_{tag}" if tag else "")
            + ".pkl")
    return base / stem


# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

def run_sweep(args: argparse.Namespace) -> Path:
    native = _load_srw_native()
    if native is None:
        raise RuntimeError(
            "_srw_native extension not built; run `make` under "
            "lrgsglib/src/lrgsglib/statsys/SignedRW/ccore/native/."
        )
    signal.signal(signal.SIGINT, _handle_stop_signal)
    signal.signal(signal.SIGTERM, _handle_stop_signal)

    p_list = np.asarray(args.pflip_list, dtype=np.float64)
    n_p   = len(p_list)
    n_d   = args.n_disorder
    n_pairs = args.n_pairs

    # Distance bins (torus Chebyshev).  Logarithmic-ish bucketing from near-neighbor
    # up to the torus half-diameter, with the top bin capturing all antipodal pairs.
    r_max = args.side // 2
    r_edges = np.unique(np.clip(np.round(
        np.geomspace(1.0, max(r_max, 2.0), num=args.n_r_bins + 1)
    ).astype(np.int64), 1, r_max + 1))
    # Ensure the leftmost edge is 0 (include all pairs) and rightmost spans r_max+1.
    r_edges = np.concatenate([[0], r_edges[1:], [r_max + 1]])
    r_edges = np.unique(r_edges)
    n_r_bins = len(r_edges) - 1

    p_touch_mean = np.full((n_p, n_d), np.nan, dtype=np.float64)
    p_touch_std  = np.full((n_p, n_d), np.nan, dtype=np.float64)
    jaccard_mean = np.full((n_p, n_d), np.nan, dtype=np.float64)
    jaccard_std  = np.full((n_p, n_d), np.nan, dtype=np.float64)
    bub_cov_mean = np.full((n_p, n_d), np.nan, dtype=np.float64)
    p_touch_r    = np.full((n_p, n_d, n_r_bins), np.nan, dtype=np.float64)
    jaccard_r    = np.full((n_p, n_d, n_r_bins), np.nan, dtype=np.float64)
    counts_r     = np.zeros((n_p, n_d, n_r_bins), dtype=np.int64)
    t_seconds    = np.zeros((n_p, n_d), dtype=np.float64)
    done_mask    = np.zeros((n_p, n_d), dtype=bool)

    proto = Lattice2D(args.side, geo=args.geo, pbc=args.pbc, pflip=0.0,
                      seed=args.seed_disorder, sgpathn=args.workdir)
    out_path = _output_path(proto, args.n_walkers, n_pairs, args.tag)

    if out_path.exists():
        try:
            with out_path.open("rb") as fh:
                prev = pickle.load(fh)
            if (prev["p_touch_mean"].shape == p_touch_mean.shape
                    and np.allclose(prev["meta"]["p_list"], p_list)):
                p_touch_mean = prev["p_touch_mean"].copy()
                p_touch_std  = prev["p_touch_std"].copy()
                jaccard_mean = prev["jaccard_mean"].copy()
                jaccard_std  = prev["jaccard_std"].copy()
                bub_cov_mean = prev["bub_cov_mean"].copy()
                t_seconds    = np.asarray(prev["t_seconds"]).copy()
                done_mask    = np.asarray(prev["done_mask"]).copy()
                # Distance-binned arrays are optional (older pickles may lack them).
                prev_ptr = prev.get("p_touch_r")
                prev_jcr = prev.get("jaccard_r")
                prev_cr  = prev.get("counts_r")
                if (prev_ptr is not None and prev_ptr.shape == p_touch_r.shape):
                    p_touch_r = np.asarray(prev_ptr).copy()
                if (prev_jcr is not None and prev_jcr.shape == jaccard_r.shape):
                    jaccard_r = np.asarray(prev_jcr).copy()
                if (prev_cr is not None and prev_cr.shape == counts_r.shape):
                    counts_r = np.asarray(prev_cr).copy()
                if args.verbose and done_mask.any():
                    n_done = int(done_mask.sum())
                    print(f"[resume] {n_done}/{n_p*n_d} slots loaded "
                          f"from {out_path}", flush=True)
        except Exception as exc:
            if args.verbose:
                print(f"[resume] could not reuse {out_path}: {exc}", flush=True)

    meta = dict(
        args=vars(args),
        p_list=p_list.tolist(),
        n_walkers=args.n_walkers,
        n_pairs=n_pairs,
        n_disorder=n_d,
        coverage=args.coverage,
        x_node=args.x_node,
        min_dist=args.min_dist,
        N=proto.N,
        r_edges=r_edges.tolist(),
    )

    def _checkpoint():
        payload = dict(
            meta=meta,
            p_touch_mean=p_touch_mean, p_touch_std=p_touch_std,
            jaccard_mean=jaccard_mean, jaccard_std=jaccard_std,
            bub_cov_mean=bub_cov_mean,
            p_touch_r=p_touch_r, jaccard_r=jaccard_r,
            counts_r=counts_r,
            t_seconds=t_seconds, done_mask=done_mask,
        )
        with out_path.open("wb") as fh:
            pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)

    rng_pairs = np.random.default_rng(args.seed_pairs)

    for i, p in enumerate(p_list):
        if _INTERRUPTED["flag"]:
            break
        t_p_start = perf_counter()
        for d in range(n_d):
            if done_mask[i, d]:
                continue
            if _INTERRUPTED["flag"]:
                break
            if args.time_budget > 0 and (perf_counter() - t_p_start) > args.time_budget:
                break
            t0 = perf_counter()
            lat = Lattice2D(args.side, geo=args.geo, pbc=args.pbc, pflip=p,
                            seed=args.seed_disorder + 1000 * d,
                            sgpathn=args.workdir)
            if p > 0:
                lat.flip_random_fract_edges()
            csr = _csr_masks(lat, args.x_node)
            stop_thresh = max(1, int(args.coverage * csr.N))

            pairs = _sample_pairs(csr.N, args.side, n_pairs, args.min_dist,
                                  rng_pairs)
            starts_A = np.repeat(pairs[:, 0], args.n_walkers).astype(np.int64)
            starts_B = np.repeat(pairs[:, 1], args.n_walkers).astype(np.int64)

            bA = _bubbles_batched(native, csr, starts_A,
                                  seed=args.seed_a + d, n_walkers=args.n_walkers,
                                  stop_thresh=stop_thresh,
                                  chunk_pairs=args.chunk_pairs)
            bB = _bubbles_batched(native, csr, starts_B,
                                  seed=args.seed_b + d, n_walkers=args.n_walkers,
                                  stop_thresh=stop_thresh,
                                  chunk_pairs=args.chunk_pairs)

            bA = bA.astype(bool)
            bB = bB.astype(bool)
            inter = (bA & bB).sum(axis=1)
            union = (bA | bB).sum(axis=1)
            touched = (inter > 0).astype(np.float64)
            jacc    = np.where(union > 0, inter / np.maximum(union, 1), 0.0)

            p_touch_mean[i, d] = float(touched.mean())
            p_touch_std[i, d]  = float(touched.std(ddof=1)) if n_pairs > 1 else 0.0
            jaccard_mean[i, d] = float(jacc.mean())
            jaccard_std[i, d]  = float(jacc.std(ddof=1)) if n_pairs > 1 else 0.0
            bub_cov_mean[i, d] = float(bA.sum(axis=1).mean() / csr.N)

            # Distance-binned observables.  Distance is in lattice units
            # (torus Chebyshev); small-r = "close" pairs, large-r = antipodal.
            dists = _pair_distances(pairs, args.side)
            t_means, t_counts = _bin_by_distance(touched, dists, r_edges)
            j_means, _        = _bin_by_distance(jacc,    dists, r_edges)
            p_touch_r[i, d] = t_means
            jaccard_r[i, d] = j_means
            counts_r[i, d]  = t_counts

            t_seconds[i, d] = perf_counter() - t0
            done_mask[i, d] = True
            _checkpoint()

        if args.verbose and done_mask[i].any():
            ptm = np.nanmean(p_touch_mean[i])
            jcm = np.nanmean(jaccard_mean[i])
            bcm = np.nanmean(bub_cov_mean[i])
            dn  = int(done_mask[i].sum())
            print(f"p={p:.4f}  P_touch={ptm:.3f}  jaccard={jcm:.4f}  "
                  f"bubble_cov={bcm:.4f}  ({dn}/{n_d} d)",
                  flush=True)

    _checkpoint()
    return out_path


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    path = run_sweep(args)
    if args.verbose:
        print(f"done -> {path}")


if __name__ == "__main__":
    main()
