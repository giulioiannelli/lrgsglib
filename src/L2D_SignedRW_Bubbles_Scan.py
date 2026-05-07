"""Parameter-scan driver for the bubble-touch experiment.

Sweeps over the 3x3x3 grid (L × n_walkers × n_trials) requested for the
under-sampling diagnosis:

    L         ∈ {16, 32, 64}
    n_walkers ∈ {N, 2N, 10N}         (N = L²)
    n_trials  ∈ {10, 50, 100}

Each of the 27 configurations runs its own L2D_SignedRW_Bubbles.run_sweep
call with n_disorder=2 on a common p-grid (skipping p=0 per convention);
each produces its own pickle tagged with the parameters, checkpointed
after every p value.  Launch with `python … &` and let it run in the
background.
"""

from __future__ import annotations

import argparse
import pickle
import sys
from argparse import Namespace
from pathlib import Path
from time import perf_counter

import numpy as np

from L2D_SignedRW_Bubbles import run_sweep

from lrgsglib.config.const import PATHDATA


DEFAULT_P_SWEEP = [
    0.005, 0.010, 0.015, 0.020, 0.025, 0.030, 0.035, 0.040, 0.045, 0.050,
    0.060, 0.070, 0.080, 0.090, 0.100, 0.110, 0.120, 0.130, 0.140, 0.150,
    0.170, 0.200, 0.250,
]


def _build_args(L: int, n_walkers: int, n_trials: int,
                nw_label: str, p_list: list[float],
                workdir: str, n_disorder: int,
                verbose: bool) -> Namespace:
    """Assemble an argparse.Namespace matching L2D_SignedRW_Bubbles._parse_args."""
    tag = f"L{L}_nw{nw_label}_nt{n_trials}"
    return Namespace(
        side=L,
        pflip_list=list(p_list),
        geo="sqr",
        pbc=True,
        x_node="reflect",
        n_walkers=n_walkers,
        n_trials=n_trials,
        n_disorder=n_disorder,
        coverage=1.0,
        node_a=0,
        node_b=-1,
        seed_a=42,
        seed_b=43,
        seed_disorder=42,
        workdir=workdir,
        tag=tag,
        verbose=verbose,
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--sides", type=int, nargs="+", default=[16, 32, 64])
    p.add_argument("--n-walker-mults", type=int, nargs="+", default=[1, 2, 10],
                   help="Multipliers M so that n_walkers = M * N (N = L**2).")
    p.add_argument("--n-trials-grid", type=int, nargs="+", default=[10, 50, 100])
    p.add_argument("--n-disorder", type=int, default=2)
    p.add_argument("--workdir", default="bubble_scan",
                   help="Subpath under PATHDATA for the pickle outputs.")
    p.add_argument("--pflip-list", type=lambda s: [float(x) for x in s.split(",") if x],
                   default=DEFAULT_P_SWEEP)
    p.add_argument("--verbose", "-v", action=argparse.BooleanOptionalAction,
                   default=True)
    p.add_argument("--skip-existing", action=argparse.BooleanOptionalAction,
                   default=True,
                   help="Skip configurations whose pickle already exists with "
                        "every p value marked done (checkpointed).")
    args = p.parse_args()

    configs: list[Namespace] = []
    for L in args.sides:
        N = L * L
        for mult in args.n_walker_mults:
            for nt in args.n_trials_grid:
                configs.append(_build_args(
                    L=L,
                    n_walkers=mult * N,
                    n_trials=nt,
                    nw_label=f"{mult}N",
                    p_list=args.pflip_list,
                    workdir=args.workdir,
                    n_disorder=args.n_disorder,
                    verbose=args.verbose,
                ))

    def _expected_pickle(cfg: Namespace) -> Path:
        N = cfg.side * cfg.side
        stem = (f"bubble_touch_L={cfg.side}_nw={cfg.n_walkers}"
                f"_nt={cfg.n_trials}"
                + (f"_{cfg.tag}" if cfg.tag else "")
                + ".pkl")
        return (PATHDATA / cfg.workdir / "l2d_squared" / "srw"
                / f"N={N}" / "bubbles" / stem)

    def _is_complete(path: Path, expected_n_p: int) -> bool:
        if not path.exists():
            return False
        try:
            with path.open("rb") as fh:
                payload = pickle.load(fh)
        except Exception:
            return False
        done = np.asarray(payload.get("done_mask"))
        if done.ndim == 0 or done.shape[0] != expected_n_p:
            return False
        return bool(np.all(done))

    total = len(configs)
    total_start = perf_counter()
    for idx, cfg in enumerate(configs, start=1):
        banner = (
            f"\n[{idx:02d}/{total:02d}] L={cfg.side}  "
            f"n_walkers={cfg.n_walkers}  n_trials={cfg.n_trials}  "
            f"n_disorder={cfg.n_disorder}  tag={cfg.tag}"
        )
        print(banner, flush=True)
        expected = _expected_pickle(cfg)
        if args.skip_existing and _is_complete(expected, len(cfg.pflip_list)):
            print(f"    already complete, skipping  -> {expected}", flush=True)
            continue
        t0 = perf_counter()
        out_path = run_sweep(cfg)
        dt = perf_counter() - t0
        print(f"    done in {dt:.1f} s  -> {out_path}", flush=True)
    print(f"\nall {total} configs completed in "
          f"{(perf_counter() - total_start) / 60:.1f} min")


if __name__ == "__main__":
    sys.exit(main())
