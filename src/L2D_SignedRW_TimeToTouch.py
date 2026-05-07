"""First-passage walker-count-to-touch experiment (Framework B).

For each disorder realisation, one random pair ``(a, b)`` is fixed and
``n_trials`` independent first-touch samples are drawn.  Each sample
uses the native ``absorb_walker_first_touch`` kernel, which alternates
absorbing walkers from the two starts and returns the walker-pair
index at which the growing bubbles first intersect:

    n_touch(p, d, t) = min { k : bubble_A(k) ∩ bubble_B(k) ≠ ∅ }

Trials that do not touch within ``n_max`` walker-pairs are censored at
``n_max + 1`` — those are informative (disconnected-cluster
realisations) and must be kept for the Kaplan-Meier analysis.

Safeguards:
* Per-(p, d) resume via ``done_mask (n_p, n_d)`` (old ``(n_p,)``
  format is read as broadcast for backward compatibility).
* ``--time-budget`` wall-clock cap per p value.
* SIGINT/SIGTERM flushes partial results to the checkpoint.

Output
------
``<PATHDATA>/<workdir>/l2d_squared/srw/N=<N>/bubbles/time_to_touch_*.pkl``
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
    print(f"\n[interrupt: signal {signum} received; finishing current (p, d) "
          "and exiting after checkpoint]", flush=True)


from lrgsglib.graphs import Lattice2D
from lrgsglib.statsys.SignedRW._kernel import _kill_masks, signed_lattice_tables
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
    p.add_argument("--n-max", type=int, required=True,
                   help="Hard cap on walker-pairs per trial.  Trials that do "
                        "not touch within n_max pairs are censored at n_max+1.")
    p.add_argument("--step-cap", type=int, default=0,
                   help="Max steps per walker (0 -> auto 10*side**2).")
    p.add_argument("--n-trials", type=int, default=500)
    p.add_argument("--n-disorder", type=int, default=2)
    p.add_argument("--seed-ft", type=int, default=42,
                   help="Base seed for first-touch calls; per-trial seed is "
                        "seed_ft + d*1e6 + t.")
    p.add_argument("--seed-disorder", type=int, default=42)
    p.add_argument("--time-budget", type=float, default=600.0,
                   help="Wall-clock budget per p value, seconds; 0 disables.")
    p.add_argument("--workdir", "-wd", default="ttt")
    p.add_argument("--tag", default="")
    p.add_argument("--verbose", "-v", action="store_true")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Kernel helpers
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


# ---------------------------------------------------------------------------
# Output path
# ---------------------------------------------------------------------------

def _output_path(lat, n_max: int, n_trials: int, tag: str) -> Path:
    base = Path(lat.path_srw) / "bubbles"
    base.mkdir(parents=True, exist_ok=True)
    stem = (f"time_to_touch_L={lat.side1}_nmax={n_max}_nt={n_trials}"
            + (f"_{tag}" if tag else "") + ".pkl")
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

    side = args.side
    N = side * side
    step_cap = args.step_cap if args.step_cap > 0 else 10 * N

    p_list = np.asarray(args.pflip_list, dtype=np.float64)
    n_p, n_d, n_t = len(p_list), args.n_disorder, args.n_trials

    n_touch = np.full((n_p, n_d, n_t), -1, dtype=np.int64)
    t_seconds = np.zeros((n_p, n_d), dtype=np.float64)
    done_mask = np.zeros((n_p, n_d), dtype=bool)

    proto = Lattice2D(side, geo=args.geo, pbc=args.pbc, pflip=0.0,
                      seed=args.seed_disorder, sgpathn=args.workdir)
    out_path = _output_path(proto, args.n_max, n_t, args.tag)

    # Resume if compatible pickle exists (handle both old (n_p,) and new (n_p, n_d)
    # done_mask layouts).
    if out_path.exists():
        try:
            with out_path.open("rb") as fh:
                prev = pickle.load(fh)
            if (prev["n_touch"].shape == n_touch.shape
                    and np.allclose(prev["meta"]["p_list"], p_list)):
                n_touch = prev["n_touch"].copy()
                prev_t = np.asarray(prev["t_seconds"])
                if prev_t.shape == t_seconds.shape:
                    t_seconds = prev_t.copy()
                elif prev_t.shape == (n_p,):
                    t_seconds = np.broadcast_to(
                        prev_t[:, None], (n_p, n_d)
                    ).copy()
                prev_dm = np.asarray(prev["done_mask"])
                if prev_dm.shape == done_mask.shape:
                    done_mask = prev_dm.copy()
                elif prev_dm.shape == (n_p,):
                    done_mask = np.broadcast_to(
                        prev_dm[:, None], (n_p, n_d)
                    ).copy()
                if args.verbose:
                    print(f"[resume] {int(done_mask.sum())}/{n_p*n_d} (p,d) "
                          f"slots loaded from {out_path}", flush=True)
        except Exception as exc:
            if args.verbose:
                print(f"[resume] could not reuse {out_path}: {exc}", flush=True)

    # Per disorder: one random (a, b) pair, fixed across trials and resumes.
    rng_starts = np.random.default_rng(args.seed_disorder)
    start_pairs = np.stack([
        rng_starts.choice(N, size=2, replace=False)
        for _ in range(n_d)
    ])

    meta = dict(
        args=vars(args),
        p_list=p_list.tolist(),
        n_max=args.n_max, step_cap=step_cap,
        n_trials=n_t, n_disorder=n_d,
        x_node=args.x_node,
        start_pairs=start_pairs.tolist(),
        N=N,
    )

    def _checkpoint():
        payload = dict(meta=meta, n_touch=n_touch,
                       t_seconds=t_seconds, done_mask=done_mask)
        with out_path.open("wb") as fh:
            pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)

    for i, p in enumerate(p_list):
        if _INTERRUPTED["flag"]:
            break
        t_p_start = perf_counter()
        budget_exceeded = False
        for d in range(n_d):
            if done_mask[i, d]:
                continue
            if _INTERRUPTED["flag"]:
                break
            if args.time_budget > 0 and (perf_counter() - t_p_start) > args.time_budget:
                budget_exceeded = True
                break
            t0 = perf_counter()
            lat = Lattice2D(side, geo=args.geo, pbc=args.pbc, pflip=p,
                            seed=args.seed_disorder + 1000 * d,
                            sgpathn=args.workdir)
            if p > 0:
                lat.flip_random_fract_edges()
            csr = _csr_masks(lat, args.x_node)
            a_d, b_d = int(start_pairs[d, 0]), int(start_pairs[d, 1])

            for t in range(n_t):
                seed_t = args.seed_ft + d * 1_000_000 + t
                out = native.absorb_walker_first_touch(
                    csr.N, csr.neigh_indices, csr.neigh_ptr,
                    csr.kill_flat, csr.refl_flat,
                    int(a_d), int(b_d),
                    int(args.n_max), int(step_cap),
                    int(seed_t), False,
                )
                n_touch[i, d, t] = int(out[0])

            t_seconds[i, d] = perf_counter() - t0
            done_mask[i, d] = True
            _checkpoint()

        if args.verbose:
            flat = n_touch[i][done_mask[i]].ravel() if done_mask[i].any() else np.array([], dtype=np.int64)
            if flat.size:
                censored = int((flat == args.n_max + 1).sum())
                median = int(np.median(flat))
                q90 = int(np.quantile(flat, 0.90))
                done_slots = int(done_mask[i].sum())
                tag = " [budget]" if budget_exceeded else ""
                print(f"p={p:.4f}  median={median:7d}  q90={q90:7d}  "
                      f"censored={censored}/{flat.size}  ({done_slots}/{n_d} d){tag}",
                      flush=True)
            else:
                print(f"p={p:.4f}  [budget exhausted before any d slot]", flush=True)
    _checkpoint()
    return out_path


def main(argv: list[str] | None = None) -> None:
    args = _parse_args(argv)
    path = run_sweep(args)
    if args.verbose:
        print(f"done -> {path}")


if __name__ == "__main__":
    main()
