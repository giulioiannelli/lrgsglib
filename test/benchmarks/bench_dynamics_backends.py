"""
Benchmark dynamics backend comparison: Python vs pybind11 vs CuPy vs C.

Compares wall-clock time for Ising Metropolis on a 2D square lattice
across all available backends at several system sizes.

Output: results saved to test/.tmp/bench_dynamics_backends_results.pkl
"""

import pickle
import time
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

TEST_TMP_DIR = Path(__file__).parent.parent / ".tmp"
TEST_TMP_DIR.mkdir(exist_ok=True)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
NUM_TRIALS = 3
T = 2.0
STEPS = 200  # Enough for timing without being too slow for Python
SEED = 42

BACKENDS = [
    ("py_met", "Python"),
    ("pb_met", "Pybind11"),
    ("cu_met", "CuPy GPU"),
]

# System sizes: Python gets slow, so limit it
CONFIGS = [
    {"side": 8, "N": 64, "py_ok": True},
    {"side": 16, "N": 256, "py_ok": True},
    {"side": 32, "N": 1024, "py_ok": True},
    {"side": 64, "N": 4096, "py_ok": False},   # skip Python
    {"side": 128, "N": 16384, "py_ok": False},  # skip Python
]


# ---------------------------------------------------------------------------
# Benchmark function
# ---------------------------------------------------------------------------
def bench_ising(side: int, runlang: str, engine: str = "nx",
                steps: int = STEPS, seed: int = SEED) -> dict:
    """Benchmark a single Ising run."""
    from lrgsglib.statsys.IsingDynamics import IsingDynamics

    if engine == "nx":
        from lrgsglib.graphs.nx.Lattice2DNX import Lattice2DNX
        lat = Lattice2DNX(side1=side, geo="sqr", pflip=0.15, seed=seed)
    else:
        from lrgsglib.graphs.gt.Lattice2DGT import Lattice2DGT
        lat = Lattice2DGT(side1=side, geo="sqr", pflip=0.15, seed=seed)
    lat.flip_random_fract_edges()

    times = []
    for trial in range(NUM_TRIALS):
        ising = IsingDynamics(
            sg=lat, T=T, steps=steps,
            runlang=runlang, seed=seed + trial,
        )
        ising.init_ising_dynamics()

        t0 = time.perf_counter()
        ising.run(verbose=False, tqdm_on=False)
        t1 = time.perf_counter()
        times.append(t1 - t0)

    return {
        "runlang": runlang, "engine": engine, "side": side,
        "N": lat.N, "steps": steps, "times": times,
        "mean": np.mean(times), "std": np.std(times),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    results = []

    print("=" * 70)
    print("BENCHMARK: Ising Dynamics Backend Comparison")
    print(f"T={T}, steps={STEPS}, {NUM_TRIALS} trials, pflip=0.15")
    print("=" * 70)

    # Check CuPy availability
    try:
        import cupy  # noqa: F401
        has_cupy = True
        print("CuPy: available")
    except ImportError:
        has_cupy = False
        print("CuPy: not available (skipping GPU benchmarks)")

    # --- NX engine benchmarks ---
    print("\n--- NX Engine ---")
    header = f"{'Backend':12s} {'side':>5s} {'N':>6s} {'time (s)':>12s} {'per sweep':>12s}"
    print(header)
    print("-" * len(header))

    for cfg in CONFIGS:
        for runlang, label in BACKENDS:
            if runlang == "py_met" and not cfg["py_ok"]:
                continue
            if runlang == "cu_met" and not has_cupy:
                continue

            try:
                r = bench_ising(cfg["side"], runlang, engine="nx")
                results.append(r)
                per_sweep = r["mean"] / STEPS
                print(f"{label:12s} {cfg['side']:5d} {r['N']:6d} "
                      f"{r['mean']:10.4f}s  {per_sweep:.6f}s")
            except Exception as e:
                print(f"{label:12s} {cfg['side']:5d} {cfg['N']:6d}  FAILED: {e}")

    # --- GT engine benchmarks (pb and cu only) ---
    print("\n--- GT Engine (pb_met and cu_met only) ---")
    print(header)
    print("-" * len(header))

    gt_backends = [("pb_met", "Pybind11"), ("cu_met", "CuPy GPU")]

    for cfg in CONFIGS:
        for runlang, label in gt_backends:
            if runlang == "cu_met" and not has_cupy:
                continue

            try:
                r = bench_ising(cfg["side"], runlang, engine="gt")
                results.append(r)
                per_sweep = r["mean"] / STEPS
                print(f"{label:12s} {cfg['side']:5d} {r['N']:6d} "
                      f"{r['mean']:10.4f}s  {per_sweep:.6f}s")
            except Exception as e:
                print(f"{label:12s} {cfg['side']:5d} {cfg['N']:6d}  FAILED: {e}")

    # --- Speedup summary ---
    print("\n" + "=" * 70)
    print("SPEEDUP vs PYTHON (NX engine)")
    print("=" * 70)

    for cfg in CONFIGS:
        if not cfg["py_ok"]:
            continue
        nx_results = [r for r in results
                      if r["side"] == cfg["side"] and r["engine"] == "nx"]
        py_time = next((r["mean"] for r in nx_results
                        if r["runlang"] == "py_met"), None)
        if py_time is None:
            continue

        print(f"\n  N={cfg['N']} (side={cfg['side']}):")
        for r in nx_results:
            if r["runlang"] == "py_met":
                continue
            speedup = py_time / r["mean"] if r["mean"] > 0 else float("inf")
            print(f"    {r['runlang']:8s}: {speedup:6.1f}x faster")

    # --- SA and PT quick test ---
    print("\n" + "=" * 70)
    print("QUICK SA/PT BACKEND TEST (side=16)")
    print("=" * 70)

    sa_pt_backends = [
        ("pb_sa", "Pybind11 SA"),
        ("pb_pt", "Pybind11 PT"),
    ]
    if has_cupy:
        sa_pt_backends.extend([
            ("cu_sa", "CuPy SA"),
            ("cu_pt", "CuPy PT"),
        ])

    from lrgsglib.graphs.nx.Lattice2DNX import Lattice2DNX
    from lrgsglib.statsys.IsingDynamics import IsingDynamics

    lat16 = Lattice2DNX(side1=16, geo="sqr", pflip=0.15, seed=42)
    lat16.flip_random_fract_edges()

    for runlang, label in sa_pt_backends:
        try:
            kw = dict(sg=lat16, runlang=runlang, seed=42)
            if "sa" in runlang:
                kw.update(sa_enabled=True, T_init=5.0, T_final=0.1,
                          n_temperatures=30, steps_per_T=10)
            else:
                kw.update(pt_enabled=True, n_replicas=4, T_min=0.5,
                          T_max=4.0, steps_per_exchange=10, n_exchanges=20)

            ising = IsingDynamics(**kw)
            ising.init_ising_dynamics()
            t0 = time.perf_counter()
            ising.run(verbose=False, tqdm_on=False)
            t1 = time.perf_counter()
            print(f"  {label:16s}: {t1-t0:.4f}s")
            results.append({
                "runlang": runlang, "engine": "nx", "side": 16,
                "N": 256, "steps": "sa/pt", "times": [t1-t0],
                "mean": t1-t0, "std": 0.0,
            })
        except Exception as e:
            print(f"  {label:16s}: FAILED ({e})")

    # Save results
    out_file = TEST_TMP_DIR / "bench_dynamics_backends_results.pkl"
    with open(out_file, "wb") as f:
        pickle.dump(results, f)
    print(f"\nResults saved to {out_file}")


if __name__ == "__main__":
    main()
