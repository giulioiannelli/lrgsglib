"""
Benchmark graph generation speed: graph-tool (GT) vs NetworkX (NX).

Compares construction time for Lattice2D, ErdosRenyi, and BarabasiAlbert
across both engines at several system sizes.

Output: results saved to test/.tmp/bench_gt_vs_nx_results.pkl
"""

import pickle
import time
from pathlib import Path

import numpy as np

TEST_TMP_DIR = Path(__file__).parent.parent / ".tmp"
TEST_TMP_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
NUM_TRIALS = 3

LATTICE_CONFIGS = [
    {"side1": 16, "geo": "sqr", "pflip": 0.2},
    {"side1": 32, "geo": "sqr", "pflip": 0.2},
    {"side1": 64, "geo": "sqr", "pflip": 0.2},
    {"side1": 128, "geo": "sqr", "pflip": 0.2},
]

ER_CONFIGS = [
    {"n": 100, "p": 0.1, "pflip": 0.2},
    {"n": 500, "p": 0.05, "pflip": 0.2},
    {"n": 1000, "p": 0.02, "pflip": 0.2},
    {"n": 2000, "p": 0.01, "pflip": 0.2},
]

BA_CONFIGS = [
    {"n": 100, "m": 3, "pflip": 0.2},
    {"n": 500, "m": 3, "pflip": 0.2},
    {"n": 1000, "m": 3, "pflip": 0.2},
    {"n": 2000, "m": 3, "pflip": 0.2},
]


# ---------------------------------------------------------------------------
# Benchmarking functions
# ---------------------------------------------------------------------------
def bench_lattice(engine: str, cfg: dict, seed: int = 42) -> dict:
    """Benchmark Lattice2D creation."""
    if engine == "nx":
        from lrgsglib.graphs.nx.Lattice2DNX import Lattice2DNX as Cls
    else:
        from lrgsglib.graphs.gt.Lattice2DGT import Lattice2DGT as Cls

    times = []
    N = None
    for trial in range(NUM_TRIALS):
        t0 = time.perf_counter()
        g = Cls(side1=cfg["side1"], geo=cfg["geo"], pflip=cfg["pflip"],
                seed=seed + trial)
        g.flip_random_fract_edges()
        t1 = time.perf_counter()
        times.append(t1 - t0)
        N = g.N

    return {
        "engine": engine, "graph": "Lattice2D", "N": N,
        "config": cfg, "times": times,
        "mean": np.mean(times), "std": np.std(times),
    }


def bench_erdos_renyi(engine: str, cfg: dict, seed: int = 42) -> dict:
    """Benchmark ErdosRenyi creation."""
    if engine == "nx":
        from lrgsglib.graphs.nx.ErdosRenyiNX import ErdosRenyiNX as Cls
    else:
        from lrgsglib.graphs.gt.ErdosRenyiGT import ErdosRenyiGT as Cls

    times = []
    N = None
    for trial in range(NUM_TRIALS):
        t0 = time.perf_counter()
        g = Cls(n=cfg["n"], p=cfg["p"], pflip=cfg["pflip"],
                seed=seed + trial)
        g.flip_random_fract_edges()
        t1 = time.perf_counter()
        times.append(t1 - t0)
        N = g.N

    return {
        "engine": engine, "graph": "ErdosRenyi", "N": N,
        "config": cfg, "times": times,
        "mean": np.mean(times), "std": np.std(times),
    }


def bench_barabasi_albert(engine: str, cfg: dict, seed: int = 42) -> dict:
    """Benchmark BarabasiAlbert creation."""
    if engine == "nx":
        from lrgsglib.graphs.nx.BarabasiAlbertNX import BarabasiAlbertNX as Cls
    else:
        from lrgsglib.graphs.gt.BarabasiAlbertGT import BarabasiAlbertGT as Cls

    times = []
    N = None
    for trial in range(NUM_TRIALS):
        t0 = time.perf_counter()
        g = Cls(n=cfg["n"], m=cfg["m"], pflip=cfg["pflip"],
                seed=seed + trial)
        g.flip_random_fract_edges()
        t1 = time.perf_counter()
        times.append(t1 - t0)
        N = g.N

    return {
        "engine": engine, "graph": "BarabasiAlbert", "N": N,
        "config": cfg, "times": times,
        "mean": np.mean(times), "std": np.std(times),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    import warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    results = []
    engines = ["nx", "gt"]

    print("=" * 70)
    print("BENCHMARK: Graph Generation — GT vs NX")
    print("=" * 70)

    # --- Lattice2D ---
    print("\n--- Lattice2D ---")
    for cfg in LATTICE_CONFIGS:
        for eng in engines:
            try:
                r = bench_lattice(eng, cfg)
                results.append(r)
                print(f"  {eng:3s} side={cfg['side1']:4d}  N={r['N']:6d}  "
                      f"{r['mean']:.4f} +/- {r['std']:.4f} s")
            except Exception as e:
                print(f"  {eng:3s} side={cfg['side1']:4d}  FAILED: {e}")

    # --- ErdosRenyi ---
    print("\n--- ErdosRenyi ---")
    for cfg in ER_CONFIGS:
        for eng in engines:
            try:
                r = bench_erdos_renyi(eng, cfg)
                results.append(r)
                print(f"  {eng:3s} n={cfg['n']:5d}  N={r['N']:6d}  "
                      f"{r['mean']:.4f} +/- {r['std']:.4f} s")
            except Exception as e:
                print(f"  {eng:3s} n={cfg['n']:5d}  FAILED: {e}")

    # --- BarabasiAlbert ---
    print("\n--- BarabasiAlbert ---")
    for cfg in BA_CONFIGS:
        for eng in engines:
            try:
                r = bench_barabasi_albert(eng, cfg)
                results.append(r)
                print(f"  {eng:3s} n={cfg['n']:5d}  N={r['N']:6d}  "
                      f"{r['mean']:.4f} +/- {r['std']:.4f} s")
            except Exception as e:
                print(f"  {eng:3s} n={cfg['n']:5d}  FAILED: {e}")

    # --- Summary: speedup ratios ---
    print("\n" + "=" * 70)
    print("SPEEDUP SUMMARY (NX time / GT time)")
    print("=" * 70)
    from itertools import groupby

    for graph_type in ["Lattice2D", "ErdosRenyi", "BarabasiAlbert"]:
        typed = [r for r in results if r["graph"] == graph_type]
        # Group by N
        by_n: dict[int, dict[str, float]] = {}
        for r in typed:
            by_n.setdefault(r["N"], {})[r["engine"]] = r["mean"]
        print(f"\n  {graph_type}:")
        for n in sorted(by_n):
            vals = by_n[n]
            if "nx" in vals and "gt" in vals:
                ratio = vals["nx"] / vals["gt"] if vals["gt"] > 0 else float("inf")
                print(f"    N={n:6d}: NX={vals['nx']:.4f}s  GT={vals['gt']:.4f}s  "
                      f"speedup={ratio:.2f}x")

    # Save results
    out_file = TEST_TMP_DIR / "bench_gt_vs_nx_results.pkl"
    with open(out_file, "wb") as f:
        pickle.dump(results, f)
    print(f"\nResults saved to {out_file}")


if __name__ == "__main__":
    main()
