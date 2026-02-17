"""
Benchmark pybind11 vs C subprocess Ising backends on 3D lattices.

Compares wall-clock time for SA (simulated annealing) across:
- pb_sa (pybind11 native) vs C3B (C subprocess SA)
- NX and GT graph engines
- Lattice sizes L=8, 16, 32

Metrics: wall time, final energy.
Output: test/.tmp/bench_ising_pb_vs_c_results.pkl
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
PFLIP = 0.15
SEED = 42

SA_PARAMS = dict(
    sa_enabled=True,
    T_init=5.0,
    T_final=0.01,
    n_temperatures=50,
    steps_per_T=100,
    cooling_schedule="geometric",
    cooling_rate=0.9,
)

BACKENDS = [
    ("pb_sa", "Pybind11 SA"),
    ("C3B", "C subprocess SA"),
]

CONFIGS = [
    {"dim": 8, "label": "8^3"},
    {"dim": 16, "label": "16^3"},
    {"dim": 32, "label": "32^3"},
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _make_lattice(dim: int, engine: str, seed: int):
    """Create a 3D SC lattice with given engine."""
    if engine == "nx":
        from lrgsglib.graphs.nx.Lattice3DNX import Lattice3DNX
        lat = Lattice3DNX(dim=dim, geo="sc", pflip=PFLIP, seed=seed)
    else:
        from lrgsglib.graphs.gt.Lattice3DGT import Lattice3DGT
        lat = Lattice3DGT(dim=dim, geo="sc", pflip=PFLIP, seed=seed)
    lat.flip_random_fract_edges()
    return lat


def bench_sa(dim: int, runlang: str, engine: str,
             seed: int = SEED) -> dict:
    """Benchmark a single SA run."""
    from lrgsglib.statsys.IsingDynamics import IsingDynamics

    lat = _make_lattice(dim, engine, seed)

    times = []
    energies = []
    for trial in range(NUM_TRIALS):
        ising = IsingDynamics(
            sg=lat, runlang=runlang, seed=seed + trial, **SA_PARAMS,
        )
        ising.init_ising_dynamics()

        t0 = time.perf_counter()
        ising.run(sa_mode=True, verbose=False, tqdm_on=False)
        t1 = time.perf_counter()

        times.append(t1 - t0)
        final_e = float(ising.ene[-1]) if hasattr(ising, "ene") and len(ising.ene) else np.nan
        energies.append(final_e)

        # Clean up C files
        if hasattr(ising, "remove_run_c_files"):
            ising.remove_run_c_files()
        if hasattr(lat, "remove_exported_files"):
            lat.remove_exported_files()

    return {
        "runlang": runlang, "engine": engine, "dim": dim,
        "N": lat.N, "times": times,
        "mean_time": np.mean(times), "std_time": np.std(times),
        "energies": energies,
        "mean_energy": np.mean(energies),
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    results = []

    print("=" * 70)
    print("BENCHMARK: pybind11 vs C subprocess SA on 3D Lattices")
    print(f"pflip={PFLIP}, SA: T_init={SA_PARAMS['T_init']}, "
          f"T_final={SA_PARAMS['T_final']}, "
          f"n_T={SA_PARAMS['n_temperatures']}, "
          f"steps/T={SA_PARAMS['steps_per_T']}")
    print(f"Trials: {NUM_TRIALS}")
    print("=" * 70)

    for engine_label, engine in [("NX", "nx"), ("GT", "gt")]:
        print(f"\n--- {engine_label} Engine ---")
        header = (f"{'Backend':16s} {'dim':>5s} {'N':>6s} "
                  f"{'time (s)':>10s} {'E_final':>10s}")
        print(header)
        print("-" * len(header))

        for cfg in CONFIGS:
            for runlang, label in BACKENDS:
                try:
                    r = bench_sa(cfg["dim"], runlang, engine)
                    results.append(r)
                    print(f"{label:16s} {cfg['dim']:5d} {r['N']:6d} "
                          f"{r['mean_time']:10.4f} {r['mean_energy']:10.2f}")
                except Exception as e:
                    print(f"{label:16s} {cfg['dim']:5d} {'?':>6s} "
                          f"FAILED: {e}")

    # --- Speedup summary ---
    print("\n" + "=" * 70)
    print("SPEEDUP: C3B vs pb_sa")
    print("=" * 70)

    for engine in ("nx", "gt"):
        eng_results = [r for r in results if r["engine"] == engine]
        print(f"\n  Engine: {engine.upper()}")
        for cfg in CONFIGS:
            dim_results = [r for r in eng_results if r["dim"] == cfg["dim"]]
            pb_time = next((r["mean_time"] for r in dim_results
                            if r["runlang"] == "pb_sa"), None)
            c_time = next((r["mean_time"] for r in dim_results
                           if r["runlang"] == "C3B"), None)
            if pb_time and c_time:
                ratio = pb_time / c_time if c_time > 0 else float("inf")
                winner = "C3B" if ratio > 1 else "pb_sa"
                print(f"    {cfg['label']:>5s}: pb_sa={pb_time:.4f}s, "
                      f"C3B={c_time:.4f}s → {winner} is "
                      f"{max(ratio, 1/ratio):.1f}x faster")

    # Save results
    out_file = TEST_TMP_DIR / "bench_ising_pb_vs_c_results.pkl"
    with open(out_file, "wb") as f:
        pickle.dump(results, f)
    print(f"\nResults saved to {out_file}")


if __name__ == "__main__":
    main()
