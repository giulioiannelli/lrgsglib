#!/usr/bin/env python
"""
Benchmark MCG_SlaplSpect.py with different backends (scipy vs cupy).

Tests eigenvalue computation performance with increasing MultiplicativeCascadeGraph
sizes to quantify cupy speedup on GPU vs scipy CPU implementation.

Output files stored in lrgsglib/test/.tmp/
"""

import subprocess
import time
from pathlib import Path
import sys
import pickle
import numpy as np

# Test output directory
TEST_TMP_DIR = Path(__file__).parent / ".tmp"
TEST_TMP_DIR.mkdir(exist_ok=True)

# Test configurations: (iterations, fraction, num_averages)
# Increasing complexity from small to medium-large graphs
TEST_CONFIGS = [
    (2, 0.5, 3, "Small: ~16 nodes"),
    (3, 0.5, 3, "Medium-small: ~64 nodes"),
    (4, 0.6, 3, "Medium: ~256 nodes"),
    (5, 0.6, 3, "Medium-large: ~1024 nodes"),
    (5, 0.7, 3, "Large: ~1024+ nodes"),
]

# Fixed parameters
P1, P2, P3, P4 = 0.8, 0.6, 0.6, 0.8
PFLIP = 0.1
MODE = "eigvals"


def get_script_path():
    """Get path to MCG_SlaplSpect.py"""
    return Path("lrgsglib/src/MCG_SlaplSpect.py")


def cleanup_output_files(iterations, fraction):
    """Remove output files from previous runs."""
    base_path = Path(f"data/msg_mc/spect/i={iterations}_f={fraction}_P={P1}_{P2}_{P3}_{P4}_regular")
    if base_path.exists():
        for pkl_file in base_path.glob("mc_eigvals_p=0.1_na=3*.pkl"):
            pkl_file.unlink()
        print(f"  Cleaned up previous results for i={iterations}, f={fraction}")


def run_mcg_slaplspect(backend, iterations, fraction, num_averages):
    """
    Run MCG_SlaplSpect.py with specified backend and parameters.

    Returns:
        tuple: (execution_time_seconds, avg_num_nodes)
    """
    script_path = get_script_path()

    cmd = [
        "python", str(script_path),
        str(P1), str(P2), str(P3), str(P4),
        str(fraction), str(iterations),
        "--mode", MODE,
        "-p", str(PFLIP),
        "--number_of_averages", str(num_averages),
        "--backend", backend,
    ]

    print(f"    Running with {backend} backend...", end=" ", flush=True)

    start_time = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    end_time = time.time()

    execution_time = end_time - start_time

    # Read output file to get actual graph sizes
    base_path = Path(f"data/msg_mc/spect/i={iterations}_f={fraction}_P={P1}_{P2}_{P3}_{P4}_regular")
    metadata_file = base_path / f"mc_eigvals_p={PFLIP}_na={num_averages}_metadata.pkl"

    avg_nodes = 0
    if metadata_file.exists():
        with open(metadata_file, 'rb') as f:
            metadata = pickle.load(f)
            avg_nodes = np.mean([m['num_nodes'] for m in metadata])

    print(f"Done in {execution_time:.3f}s (avg {avg_nodes:.1f} nodes)")

    return execution_time, avg_nodes


def run_benchmark():
    """Run full benchmark suite."""
    print("=" * 80)
    print("MCG_SlaplSpect Backend Benchmark: scipy vs cupy")
    print("=" * 80)
    print(f"\nFixed parameters:")
    print(f"  Probability matrix: P = [{P1}, {P2}, {P3}, {P4}]")
    print(f"  Edge flip probability: pflip = {PFLIP}")
    print(f"  Mode: {MODE}")
    print()

    results = []

    for iterations, fraction, num_averages, description in TEST_CONFIGS:
        print(f"\n{description}")
        print(f"  Parameters: iterations={iterations}, fraction={fraction}, num_averages={num_averages}")

        # Clean up previous runs
        cleanup_output_files(iterations, fraction)

        # Test scipy backend
        try:
            scipy_time, avg_nodes = run_mcg_slaplspect(
                "scipy", iterations, fraction, num_averages
            )
        except subprocess.CalledProcessError as e:
            print(f"    ERROR: scipy backend failed")
            print(f"    {e.stderr}")
            scipy_time = None
            avg_nodes = 0

        # Clean up for cupy test
        cleanup_output_files(iterations, fraction)

        # Test cupy backend
        try:
            cupy_time, _ = run_mcg_slaplspect(
                "cupy", iterations, fraction, num_averages
            )
        except subprocess.CalledProcessError as e:
            print(f"    ERROR: cupy backend failed")
            print(f"    {e.stderr}")
            cupy_time = None
        except Exception as e:
            print(f"    WARNING: cupy backend unavailable ({e})")
            cupy_time = None

        # Calculate speedup
        if scipy_time and cupy_time:
            speedup = scipy_time / cupy_time
        else:
            speedup = None

        results.append({
            'iterations': iterations,
            'fraction': fraction,
            'num_averages': num_averages,
            'description': description,
            'avg_nodes': avg_nodes,
            'scipy_time': scipy_time,
            'cupy_time': cupy_time,
            'speedup': speedup
        })

    return results


def print_results_table(results):
    """Print formatted results table."""
    print("\n" + "=" * 80)
    print("BENCHMARK RESULTS")
    print("=" * 80)
    print()
    print(f"{'Config':<25} {'Avg Nodes':<12} {'Scipy (s)':<12} {'Cupy (s)':<12} {'Speedup':<10}")
    print("-" * 80)

    for r in results:
        config_str = f"i={r['iterations']}, f={r['fraction']}"

        scipy_str = f"{r['scipy_time']:.3f}" if r['scipy_time'] else "FAILED"
        cupy_str = f"{r['cupy_time']:.3f}" if r['cupy_time'] else "N/A"

        if r['speedup']:
            speedup_str = f"{r['speedup']:.2f}x"
            if r['speedup'] > 1:
                speedup_str += " ⚡"
        else:
            speedup_str = "-"

        print(f"{config_str:<25} {r['avg_nodes']:<12.1f} {scipy_str:<12} {cupy_str:<12} {speedup_str:<10}")

    print()

    # Summary statistics
    valid_speedups = [r['speedup'] for r in results if r['speedup']]
    if valid_speedups:
        print("SUMMARY:")
        print(f"  Average speedup: {np.mean(valid_speedups):.2f}x")
        print(f"  Max speedup: {np.max(valid_speedups):.2f}x")
        print(f"  Min speedup: {np.min(valid_speedups):.2f}x")
        print()

        # Check if speedup increases with size
        if len(valid_speedups) >= 2:
            increasing = all(valid_speedups[i] <= valid_speedups[i+1]
                           for i in range(len(valid_speedups)-1))
            if increasing:
                print("  ✓ Speedup increases with graph size (as expected)")
            else:
                print("  ⚠ Speedup does not consistently increase with size")
    else:
        print("WARNING: No valid speedup measurements obtained.")
        print("This usually means cupy is not available or failed to run.")

    print()


def check_cupy_available():
    """Check if cupy is installed and usable."""
    try:
        import cupy
        print(f"✓ CuPy version {cupy.__version__} detected")
        # Try a simple operation
        x = cupy.array([1, 2, 3])
        print(f"✓ CuPy GPU operations working")
        return True
    except ImportError:
        print("✗ CuPy not installed")
        print("  Install with: conda install -c conda-forge cupy")
        return False
    except Exception as e:
        print(f"✗ CuPy installed but not working: {e}")
        return False


def main():
    """Main entry point."""
    print("Checking CuPy availability...")
    cupy_available = check_cupy_available()
    print()

    if not cupy_available:
        print("WARNING: Proceeding anyway to test scipy backend.")
        print("Cupy comparisons will be skipped.\n")

    # Run benchmarks
    results = run_benchmark()

    # Print results
    print_results_table(results)

    # Save results
    output_file = TEST_TMP_DIR / "bench_mcg_slaplspect_backends_results.pkl"
    with open(output_file, 'wb') as f:
        pickle.dump(results, f)
    print(f"Results saved to: {output_file}")


if __name__ == "__main__":
    main()
