"""Bubble-touch communication-probability sweep on a 2D signed lattice.

Each trial runs `n_walkers` absorbing-walker simulations from node `a` and
another `n_walkers` from node `b`; the **bubbles** are the unions of their
visited sites.  A trial succeeds iff the two bubbles intersect.  For each
disorder realisation and each `p` in the sweep we compute:

* ``P_touch(p, d) = (# successes) / n_trials``
* ``bubble_cov(p, d) = mean fraction of the lattice covered by bubble A``

Averaged over disorder realisations these give the communication-probability
order parameter `⟨P_touch⟩` vs `p`.

The program mirrors the run_absorb_walker pybind11 backend directly; the
expensive inner loop happens in C.  Designed to be run in the background
(``python L2D_SignedRW_Bubbles.py … &``); partial results are checkpointed
to the output pickle after every `p` value so a kill at any time leaves a
valid partial pickle.

Usage
-----
python lrgsglib/src/L2D_SignedRW_Bubbles.py \
    --side 64 --pflip-list 0.01,0.02,...,0.25 \
    --n-walkers 500 --n-trials 100 --n-disorder 2 \
    --workdir bubble_L=64

The output file lands under
``<PATHDATA>/<workdir>/l2d_squared/srw/N=<N>/bubbles/<stem>.pkl``.
"""

from __future__ import annotations

import argparse
import pickle
import signal
import sys
from dataclasses import asdict, dataclass
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
    p.add_argument("--side", "-s", type=int, required=True,
                   help="Square-lattice side length.")
    p.add_argument("--pflip-list", required=True, type=_parse_float_list,
                   help="Comma-separated list of flip probabilities p.")
    p.add_argument("--geo", default="sqr")
    p.add_argument("--pbc", action=argparse.BooleanOptionalAction, default=True)
    p.add_argument("--x-node", choices=["reflect", "absorb"], default="reflect",
                   help="Unfrustrated-X-node behaviour for absorbing walker.")
    p.add_argument("--n-walkers", "-nw", type=int, default=500,
                   help="Walkers per bubble (ensemble) per trial.")
    p.add_argument("--n-trials", "-nt", type=int, default=100,
                   help="Independent trials per disorder realisation.")
    p.add_argument("--n-disorder", "-nd", type=int, default=2,
                   help="Disorder realisations averaged.")
    p.add_argument("--coverage", type=float, default=1.0,
                   help="Coverage-stop fraction (1.0 = walkers run until death).")
    p.add_argument("--node-a", type=int, default=0,
                   help="Absolute node index for bubble A (default: node 0).")
    p.add_argument("--node-b", type=int, default=-1,
                   help="Absolute node index for bubble B "
                        "(default: lattice centre (side//2)*(side+1)).")
    p.add_argument("--seed-a", type=int, default=42)
    p.add_argument("--seed-b", type=int, default=43)
    p.add_argument("--seed-disorder", type=int, default=42,
                   help="Base seed for disorder-realisation sign flips; "
                        "per-realisation seed is seed_disorder + 1000*d.")
    p.add_argument("--time-budget", type=float, default=0.0,
                   help="Wall-clock budget per p value, seconds (0=disabled).")
    p.add_argument("--workdir", "-wd", default="",
                   help="Subpath under PATHDATA (e.g. 'bubble_L=64').")
    p.add_argument("--tag", default="",
                   help="Extra suffix on the output filename stem.")
    p.add_argument("--verbose", "-v", action="store_true")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Native-kernel helpers
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


def _bubbles_from_site(native, csr: _CSRMasks, *,
                       node: int, seed: int, n_walkers: int,
                       n_trials: int, stop_thresh: int) -> np.ndarray:
    total = n_walkers * n_trials
    starts = np.full(total, int(node), dtype=np.int64)
    out = native.absorb_walker_sampling(
        csr.N,
        csr.neigh_indices, csr.neigh_ptr,
        csr.kill_flat, csr.refl_flat,
        starts,
        int(stop_thresh),
        int(seed),
        int(n_walkers),   # walkers_per_trial -> kernel emits (n_trials, N) bubble bitmap
        False,             # track_unique=False: skip per-walker bitmap (fast path)
    )
    return out[5]  # (n_trials, N) uint8


# ---------------------------------------------------------------------------
# Output path
# ---------------------------------------------------------------------------

def _output_path(lat, n_walkers: int, n_trials: int, tag: str) -> Path:
    base = Path(lat.path_srw) / "bubbles"
    base.mkdir(parents=True, exist_ok=True)
    stem = (f"bubble_touch_L={lat.side1}_nw={n_walkers}_nt={n_trials}"
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
    n_p = len(p_list)
    n_d = args.n_disorder
    #
    p_touch   = np.full((n_p, n_d), np.nan, dtype=np.float64)
    bub_cov   = np.full((n_p, n_d), np.nan, dtype=np.float64)
    t_seconds = np.zeros((n_p, n_d), dtype=np.float64)
    done_mask = np.zeros((n_p, n_d), dtype=bool)
    #
    # Build one lattice to discover path_srw; the per-p loop rebuilds with
    # the right pflip.
    proto = Lattice2D(args.side, geo=args.geo, pbc=args.pbc, pflip=0.0,
                      seed=args.seed_disorder, sgpathn=args.workdir)
    node_a = args.node_a
    node_b = (args.node_b if args.node_b >= 0
              else (args.side // 2) * (args.side + 1))
    out_path = _output_path(proto, args.n_walkers, args.n_trials, args.tag)

    # Resume from an existing checkpoint if one matches the current shape.
    if out_path.exists():
        try:
            with out_path.open("rb") as fh:
                prev = pickle.load(fh)
            if (prev["p_touch"].shape == p_touch.shape
                    and np.allclose(prev["meta"]["p_list"], p_list)):
                p_touch   = prev["p_touch"].copy()
                bub_cov   = prev["bub_cov"].copy()
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
                n_resumed = int(done_mask.sum())
                if args.verbose and n_resumed:
                    print(f"[resume] {n_resumed}/{n_p*n_d} (p,d) slots "
                          f"loaded from {out_path}", flush=True)
        except Exception as exc:
            if args.verbose:
                print(f"[resume] could not reuse {out_path}: {exc}", flush=True)
    #
    meta = dict(
        args=vars(args),
        p_list=p_list.tolist(),
        n_walkers=args.n_walkers,
        n_trials=args.n_trials,
        n_disorder=n_d,
        coverage=args.coverage,
        x_node=args.x_node,
        node_a=node_a, node_b=node_b,
        N=proto.N,
    )
    #
    def _checkpoint():
        payload = dict(
            meta=meta,
            p_touch=p_touch, bub_cov=bub_cov,
            t_seconds=t_seconds, done_mask=done_mask,
        )
        with out_path.open("wb") as fh:
            pickle.dump(payload, fh, protocol=pickle.HIGHEST_PROTOCOL)
    #
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
            lat = Lattice2D(args.side, geo=args.geo, pbc=args.pbc, pflip=p,
                            seed=args.seed_disorder + 1000 * d,
                            sgpathn=args.workdir)
            if p > 0:
                lat.flip_random_fract_edges()
            csr = _csr_masks(lat, args.x_node)
            stop_thresh = max(1, int(args.coverage * csr.N))
            bA = _bubbles_from_site(
                native, csr,
                node=node_a, seed=args.seed_a + d,
                n_walkers=args.n_walkers, n_trials=args.n_trials,
                stop_thresh=stop_thresh,
            )
            bB = _bubbles_from_site(
                native, csr,
                node=node_b, seed=args.seed_b + d,
                n_walkers=args.n_walkers, n_trials=args.n_trials,
                stop_thresh=stop_thresh,
            )
            touched = (bA.astype(bool) & bB.astype(bool)).any(axis=1)
            p_touch[i, d] = float(touched.mean())
            bub_cov[i, d] = float(bA.sum(axis=1).mean() / csr.N)
            t_seconds[i, d] = perf_counter() - t0
            done_mask[i, d] = True
            _checkpoint()
        if args.verbose and done_mask[i].any():
            done_slots = int(done_mask[i].sum())
            ptm = np.nanmean(p_touch[i])
            tag = " [budget]" if budget_exceeded else ""
            print(f"p={p:.4f}  P_touch={ptm:.3f}  "
                  f"bubble_cov={np.nanmean(bub_cov[i]):.4f}  "
                  f"({done_slots}/{n_d} d){tag}",
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
