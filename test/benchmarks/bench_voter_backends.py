"""
Benchmark VoterModel backend / rule / sampler comparison.

Three studies, all on signed 2D lattices (and Erdos-Renyi) via the
``lrgsglib.graphs`` factories:

1. **Backend scaling** -- synchronous *linear* voter across all backends
   (``py``, ``C0``, ``pb_voter``, ``np_voter``, ``cu_voter``) at several system
   sizes; speedup vs pure Python. This is the apples-to-apples comparison
   because every backend implements the synchronous linear voter.
2. **Sampler study** -- rejection-free ``gillespie`` CTMC vs ``asynchronous``
   random-sequential, with ``absorbing_check=True`` on a *balanced* (unsigned)
   lattice where consensus is reachable. Gillespie's advantage grows near
   consensus (the low-activity coarsening tail).
3. **Rule study** -- per-rule overhead (``linear`` / ``majority`` / ``qvoter`` /
   ``nonlinear``) at fixed size/backend.

The module exposes ``run_backend_scaling``, ``run_sampler_study`` and
``run_rule_study`` so ``ipy/devtool/VoterModel_benchmarks.ipynb`` can drive the
exact same tested code. Running the file as a script executes all three and
pickles the results to ``test/.tmp/bench_voter_backends_results.pkl``.
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
# Configuration (named constants -- no inline magic numbers)
# ---------------------------------------------------------------------------
NUM_TRIALS = 3
SEED = 42

# Signed substrate: fraction of edges flipped for the (frustrated) scaling and
# rule studies; the sampler study uses a balanced graph (pflip = 0) so an
# absorbing consensus exists.
PFLIP_SIGNED = 0.15
PFLIP_BALANCED = 0.0

# Erdos-Renyi target mean degree (maps side -> N = side**2, p = k / (N - 1)).
ER_AVG_DEGREE = 6.0

# --- Study 1: backend scaling (synchronous linear voter, fixed sweeps) ---
SCALING_STEPS = 200
SCALING_CONFIGS = [
    {"side": 16, "py_ok": True},
    {"side": 32, "py_ok": True},
    {"side": 64, "py_ok": False},    # Python too slow above here
    {"side": 128, "py_ok": False},
    {"side": 256, "py_ok": False},   # GPU crossover territory
]
# Backends that all implement the synchronous linear voter.
SCALING_BACKENDS = [
    ("py", "Python", True),          # (runlang, label, needs_py_ok)
    ("C0", "C subproc", False),
    ("pb_voter", "pybind11", False),
    ("np_voter", "NumPy vec", False),
    ("cu_voter", "CuPy GPU", False),
]

# --- Study 2: sampler study (gillespie vs asynchronous, absorbing) ---
SAMPLER_STEPS = 4000            # generous cap; absorbing_check stops early
SAMPLER_CONFIGS = [16, 32, 48, 64]
SAMPLER_BACKENDS = [
    ("py", "Python"),
    ("pb_voter", "pybind11"),
    ("C0", "C subproc"),
]

# --- Study 3: rule overhead (asynchronous, fixed size/backend) ---
RULE_STUDY_SIDE = 48
RULE_STUDY_STEPS = 300
RULE_STUDY_BACKEND = "pb_voter"
RULE_STUDY_RULES = ["linear", "majority", "qvoter", "nonlinear"]

# --- Study 4: cluster-tracking overhead (gillespie, balanced lattice) ---
# Cost of track_clusters=True vs False on the single-spin gillespie path.
CLUSTER_STUDY_BACKENDS = ["py", "pb_voter"]
CLUSTER_STUDY_SIDES = [24, 32, 48]
CLUSTER_STUDY_STEPS = 1500

# Graph families exercised by the scaling sweep. C0 is NX-only, hence NX
# everywhere; lattice geometries keep contiguous 0..N-1 node labels (required by
# the index-based backends). 'er' is supported by build_graph but only when the
# instance is fully connected -- the giant-component subgraph otherwise has
# non-contiguous labels the backends reject -- so it is not in the default sweep;
# the "signed structure" axis is covered by pflip in the sampler/cluster studies.
GRAPH_FAMILIES = ["sqr", "hex"]


# ---------------------------------------------------------------------------
# Graph factory
# ---------------------------------------------------------------------------
def build_graph(family: str, side: int, pflip: float, seed: int):
    """Build a signed graph for the requested family (NX engine).

    ``side`` sets the linear lattice size, or maps to ``N = side**2`` for the
    Erdos-Renyi family. Sign flips are applied in-place.
    """
    if family in ("sqr", "tri", "hex"):
        from lrgsglib.graphs.nx import Lattice2DNX
        g = Lattice2DNX(side1=side, geo=family, pflip=pflip, seed=seed)
    elif family == "er":
        from lrgsglib.graphs.nx import ErdosRenyiNX
        n = side * side
        p = ER_AVG_DEGREE / (n - 1)
        g = ErdosRenyiNX(n=n, p=p, pflip=pflip, seed=seed)
    else:
        raise ValueError(f"unknown graph family: {family!r}")
    g.flip_random_fract_edges()
    return g


def cupy_available() -> bool:
    try:
        import cupy  # noqa: F401
        return True
    except Exception:
        return False


# ---------------------------------------------------------------------------
# Single-run timer
# ---------------------------------------------------------------------------
def bench_voter(
    family: str,
    side: int,
    runlang: str,
    *,
    rule: str = "linear",
    upd_mode: str = "synchronous",
    steps: int = SCALING_STEPS,
    pflip: float = PFLIP_SIGNED,
    absorbing_check: bool = False,
    track_clusters: bool = False,
    n_trials: int = NUM_TRIALS,
    seed: int = SEED,
) -> dict:
    """Time a single VoterModel configuration over ``n_trials`` runs.

    Returns a record with per-trial wall times, the mean/std, the actual node
    count, and (when ``absorbing_check``) the mean absorption sweep.
    """
    from lrgsglib.statsys.VoterModel import VoterModel

    times: list[float] = []
    absorbed: list[int] = []
    N = None
    for trial in range(n_trials):
        g = build_graph(family, side, pflip, seed + trial)
        N = g.N
        vm = VoterModel(
            sg=g, steps=steps, runlang=runlang,
            rule=rule, upd_mode=upd_mode,
            absorbing_check=absorbing_check,
            track_clusters=track_clusters,
            seed=seed + trial,
        )
        vm.init_voter_dynamics()
        t0 = time.perf_counter()
        vm.run(verbose=False, tqdm_on=False)
        t1 = time.perf_counter()
        times.append(t1 - t0)
        if absorbing_check and vm.absorbed_at is not None:
            absorbed.append(int(vm.absorbed_at))

    return {
        "study": None,
        "family": family, "side": side, "N": N,
        "runlang": runlang, "rule": rule, "upd_mode": upd_mode,
        "steps": steps, "pflip": pflip, "absorbing_check": absorbing_check,
        "track_clusters": track_clusters,
        "times": times, "mean": float(np.mean(times)),
        "std": float(np.std(times)),
        "absorbed_mean": (float(np.mean(absorbed)) if absorbed else None),
    }


# ---------------------------------------------------------------------------
# Study 1: backend scaling
# ---------------------------------------------------------------------------
def run_backend_scaling(
    families=GRAPH_FAMILIES, verbose: bool = True,
) -> list[dict]:
    """Synchronous linear voter across backends and sizes; speedup vs py."""
    has_cupy = cupy_available()
    results: list[dict] = []
    for family in families:
        if verbose:
            print(f"\n--- Backend scaling: family={family} ---")
            header = (f"{'Backend':12s} {'side':>5s} {'N':>7s} "
                      f"{'time (s)':>12s} {'per sweep':>12s} {'speedup':>9s}")
            print(header)
            print("-" * len(header))
        for cfg in SCALING_CONFIGS:
            py_time = None
            for runlang, label, needs_py_ok in SCALING_BACKENDS:
                if needs_py_ok and not cfg["py_ok"]:
                    continue
                if runlang == "cu_voter" and not has_cupy:
                    continue
                try:
                    r = bench_voter(
                        family, cfg["side"], runlang,
                        rule="linear", upd_mode="synchronous",
                        steps=SCALING_STEPS, pflip=PFLIP_SIGNED,
                    )
                    r["study"] = "scaling"
                    if runlang == "py":
                        py_time = r["mean"]
                    speedup = (py_time / r["mean"]
                               if py_time and r["mean"] > 0 else None)
                    r["speedup_vs_py"] = speedup
                    results.append(r)
                    if verbose:
                        per = r["mean"] / SCALING_STEPS
                        sp = f"{speedup:7.1f}x" if speedup else "      --"
                        print(f"{label:12s} {cfg['side']:5d} {r['N']:7d} "
                              f"{r['mean']:10.4f}s  {per:.6f}s  {sp}")
                except Exception as e:
                    if verbose:
                        print(f"{label:12s} {cfg['side']:5d} {'':7s} "
                              f"FAILED: {type(e).__name__}: {e}")
    return results


# ---------------------------------------------------------------------------
# Study 2: sampler study (gillespie vs asynchronous, absorbing)
# ---------------------------------------------------------------------------
def run_sampler_study(verbose: bool = True) -> list[dict]:
    """Rejection-free gillespie vs asynchronous on a balanced lattice."""
    results: list[dict] = []
    if verbose:
        print("\n--- Sampler study: gillespie vs asynchronous "
              "(balanced sqr lattice, absorbing) ---")
        header = (f"{'side':>5s} {'N':>7s} {'async (s)':>12s} "
                  f"{'gillespie (s)':>14s} {'speedup':>9s} {'absorb@':>10s}")
        print(header)
        print("-" * len(header))
    for backend, label in SAMPLER_BACKENDS:
        if verbose:
            print(f"  backend = {label} ({backend})")
        for side in SAMPLER_CONFIGS:
            try:
                r_async = bench_voter(
                    "sqr", side, backend, rule="linear",
                    upd_mode="asynchronous", steps=SAMPLER_STEPS,
                    pflip=PFLIP_BALANCED, absorbing_check=True,
                )
                r_gil = bench_voter(
                    "sqr", side, backend, rule="linear",
                    upd_mode="gillespie", steps=SAMPLER_STEPS,
                    pflip=PFLIP_BALANCED, absorbing_check=True,
                )
                r_async["study"] = "sampler"
                r_gil["study"] = "sampler"
                speedup = (r_async["mean"] / r_gil["mean"]
                           if r_gil["mean"] > 0 else None)
                r_gil["speedup_vs_async"] = speedup
                results.extend([r_async, r_gil])
                if verbose:
                    sp = f"{speedup:7.2f}x" if speedup else "      --"
                    ab = (f"{r_gil['absorbed_mean']:.0f}"
                          if r_gil["absorbed_mean"] is not None else "--")
                    print(f"{side:5d} {r_async['N']:7d} "
                          f"{r_async['mean']:10.4f}s  {r_gil['mean']:12.4f}s  "
                          f"{sp}  {ab:>10s}")
            except Exception as e:
                if verbose:
                    print(f"{side:5d}  FAILED: {type(e).__name__}: {e}")
    return results


# ---------------------------------------------------------------------------
# Study 3: rule overhead
# ---------------------------------------------------------------------------
def run_rule_study(verbose: bool = True) -> list[dict]:
    """Per-rule overhead at fixed size/backend (asynchronous schedule)."""
    results: list[dict] = []
    if verbose:
        print(f"\n--- Rule study: {RULE_STUDY_BACKEND}, side={RULE_STUDY_SIDE}, "
              f"asynchronous ---")
        header = f"{'rule':12s} {'N':>7s} {'time (s)':>12s} {'vs linear':>10s}"
        print(header)
        print("-" * len(header))
    lin_time = None
    for rule in RULE_STUDY_RULES:
        try:
            r = bench_voter(
                "sqr", RULE_STUDY_SIDE, RULE_STUDY_BACKEND,
                rule=rule, upd_mode="asynchronous",
                steps=RULE_STUDY_STEPS, pflip=PFLIP_SIGNED,
            )
            r["study"] = "rule"
            if rule == "linear":
                lin_time = r["mean"]
            ratio = (r["mean"] / lin_time if lin_time and lin_time > 0 else None)
            r["ratio_vs_linear"] = ratio
            results.append(r)
            if verbose:
                rt = f"{ratio:8.2f}x" if ratio else "      --"
                print(f"{rule:12s} {r['N']:7d} {r['mean']:10.4f}s  {rt}")
        except Exception as e:
            if verbose:
                print(f"{rule:12s}  FAILED: {type(e).__name__}: {e}")
    return results


# ---------------------------------------------------------------------------
# Study 4: cluster-tracking overhead
# ---------------------------------------------------------------------------
def run_cluster_overhead(verbose: bool = True) -> list[dict]:
    """Overhead of track_clusters=True vs False on the gillespie path."""
    results: list[dict] = []
    if verbose:
        print("\n--- Cluster-tracking overhead (gillespie, balanced sqr) ---")
        header = (f"{'Backend':10s} {'side':>5s} {'N':>7s} {'off (s)':>11s} "
                  f"{'on (s)':>11s} {'overhead':>10s}")
        print(header)
        print("-" * len(header))
    for backend in CLUSTER_STUDY_BACKENDS:
        for side in CLUSTER_STUDY_SIDES:
            try:
                off = bench_voter(
                    "sqr", side, backend, rule="linear", upd_mode="gillespie",
                    steps=CLUSTER_STUDY_STEPS, pflip=PFLIP_BALANCED,
                    absorbing_check=True, track_clusters=False,
                )
                on = bench_voter(
                    "sqr", side, backend, rule="linear", upd_mode="gillespie",
                    steps=CLUSTER_STUDY_STEPS, pflip=PFLIP_BALANCED,
                    absorbing_check=True, track_clusters=True,
                )
                off["study"] = on["study"] = "cluster"
                ratio = on["mean"] / off["mean"] if off["mean"] > 0 else None
                on["overhead_vs_untracked"] = ratio
                results.extend([off, on])
                if verbose:
                    ov = f"{ratio:8.2f}x" if ratio else "      --"
                    print(f"{backend:10s} {side:5d} {off['N']:7d} "
                          f"{off['mean']:9.4f}s  {on['mean']:9.4f}s  {ov}")
            except Exception as e:
                if verbose:
                    print(f"{backend:10s} {side:5d}  FAILED: "
                          f"{type(e).__name__}: {e}")
    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 70)
    print("BENCHMARK: VoterModel backends / rules / samplers / clusters")
    print(f"{NUM_TRIALS} trials, seed={SEED}, pflip(signed)={PFLIP_SIGNED}")
    print(f"CuPy: {'available' if cupy_available() else 'NOT available'}")
    print("=" * 70)

    results: list[dict] = []
    results += run_backend_scaling()
    results += run_sampler_study()
    results += run_rule_study()
    results += run_cluster_overhead()

    out_file = TEST_TMP_DIR / "bench_voter_backends_results.pkl"
    with open(out_file, "wb") as f:
        pickle.dump(results, f)
    print(f"\nResults saved to {out_file}")


if __name__ == "__main__":
    main()
