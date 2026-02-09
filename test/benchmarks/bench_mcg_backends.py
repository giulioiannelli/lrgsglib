"""
Benchmark script for MCG eigenvalue computation backends.

Compares performance of scipy, numpy, and cupy backends for
eigenvalue computation on Multiplicative Cascade Graphs.

System sizes: nodes = fraction × 4^iterations
"""

import time
import sys
from pathlib import Path
import numpy as np

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "lrgsglib" / "src"))

from lrgsglib.graphs.nx import MultiplicativeCascadeGraphNX as MultiplicativeCascadeGraph


# Test configurations (user specified: iter 3-5, fraction 0.7)
TEST_CONFIGS = [
    # (iterations, fraction, expected_nodes, description)
    (3, 0.7, 44, "Very small - sanity check"),
    (4, 0.7, 179, "Tiny - verify all backends work"),
    (5, 0.7, 716, "Small - baseline timing"),
]

# Extended configs (recommended for seeing cupy speedup)
EXTENDED_CONFIGS = [
    (6, 0.7, 2867, "Small but shows cupy benefits"),
    (7, 0.7, 11468, "Medium - clear cupy speedup"),
    (8, 0.5, 32768, "Large - significant cupy advantage"),
]


def check_cupy_available():
    """Check if cupy is available."""
    try:
        import cupy
        return True
    except ImportError:
        return False


def benchmark_backend(iterations, fraction, backend, num_trials=3):
    """
    Benchmark a single backend configuration.

    Parameters
    ----------
    iterations : int
        MCG iterations
    fraction : float
        Node sampling fraction
    backend : str
        Backend to use ('scipy', 'numpy', or 'cupy')
    num_trials : int
        Number of trials to average

    Returns
    -------
    dict
        Results with timing information
    """
    print(f"\n  Testing backend '{backend}' with {num_trials} trials...", end=" ")

    times_gen = []
    times_eig = []
    times_total = []
    actual_nodes = None

    try:
        for trial in range(num_trials):
            # Time graph generation
            t0 = time.time()
            G = MultiplicativeCascadeGraph(
                p1=0.8, p2=0.6, p3=0.6, p4=0.8,
                fraction=fraction,
                iterations=iterations,
                stochastic=False,
                periodic=False,
                variant='exp_clocks',
                pflip=0.0
            )
            t1 = time.time()
            time_gen = t1 - t0

            actual_nodes = G.N

            # Time eigenvalue computation
            t2 = time.time()
            G.compute_laplacian_spectrum_weigV(backend=backend)
            t3 = time.time()
            time_eig = t3 - t2

            time_total = t3 - t0

            times_gen.append(time_gen)
            times_eig.append(time_eig)
            times_total.append(time_total)

        # Compute statistics
        mean_gen = np.mean(times_gen)
        std_gen = np.std(times_gen)
        mean_eig = np.mean(times_eig)
        std_eig = np.std(times_eig)
        mean_total = np.mean(times_total)
        std_total = np.std(times_total)

        print(f"✓ {mean_total:.4f}s")

        return {
            'backend': backend,
            'success': True,
            'actual_nodes': actual_nodes,
            'time_gen_mean': mean_gen,
            'time_gen_std': std_gen,
            'time_eig_mean': mean_eig,
            'time_eig_std': std_eig,
            'time_total_mean': mean_total,
            'time_total_std': std_total,
            'num_trials': num_trials
        }

    except Exception as e:
        print(f"✗ FAILED: {e}")
        return {
            'backend': backend,
            'success': False,
            'error': str(e)
        }


def run_benchmark(configs, use_extended=False):
    """
    Run benchmark for all configurations.

    Parameters
    ----------
    configs : list
        List of (iterations, fraction, expected_nodes, description) tuples
    use_extended : bool
        Whether to include extended configurations
    """
    print("\n" + "="*70)
    print("MCG Backend Benchmark")
    print("="*70)

    # Check cupy availability
    cupy_available = check_cupy_available()
    backends = ['scipy', 'numpy']
    if cupy_available:
        backends.append('cupy')
        print("\n✓ cupy is available")
    else:
        print("\n⚠ cupy is NOT available (will only test scipy and numpy)")

    print(f"\nBackends to test: {backends}")
    print(f"Number of configurations: {len(configs)}")
    print(f"Trials per configuration: 3")

    # Run benchmarks
    all_results = []

    for iterations, fraction, expected_nodes, description in configs:
        print(f"\n" + "-"*70)
        print(f"Configuration: iterations={iterations}, fraction={fraction}")
        print(f"  Expected nodes: ~{expected_nodes}")
        print(f"  Description: {description}")

        config_results = {
            'iterations': iterations,
            'fraction': fraction,
            'expected_nodes': expected_nodes,
            'description': description,
            'backends': {}
        }

        for backend in backends:
            result = benchmark_backend(iterations, fraction, backend, num_trials=3)
            config_results['backends'][backend] = result

        all_results.append(config_results)

    # Print summary
    print("\n" + "="*70)
    print("BENCHMARK RESULTS SUMMARY")
    print("="*70)

    for config in all_results:
        print(f"\n{config['description']} (iter={config['iterations']}, frac={config['fraction']})")
        print(f"  Expected nodes: ~{config['expected_nodes']}")

        # Get actual nodes from first successful backend
        actual_nodes = None
        for backend_name, result in config['backends'].items():
            if result['success'] and actual_nodes is None:
                actual_nodes = result['actual_nodes']
                break

        if actual_nodes:
            print(f"  Actual nodes: {actual_nodes}")

        print(f"\n  {'Backend':<10} {'Total Time (s)':<20} {'Eigenvalue Time (s)':<25} {'Speedup':<10}")
        print(f"  {'-'*10} {'-'*20} {'-'*25} {'-'*10}")

        # Reference time (scipy)
        scipy_result = config['backends'].get('scipy', {})
        ref_time = scipy_result.get('time_total_mean', None) if scipy_result.get('success') else None

        for backend_name in backends:
            result = config['backends'][backend_name]

            if result['success']:
                total_mean = result['time_total_mean']
                total_std = result['time_total_std']
                eig_mean = result['time_eig_mean']
                eig_std = result['time_eig_std']

                total_str = f"{total_mean:.4f} ± {total_std:.4f}"
                eig_str = f"{eig_mean:.4f} ± {eig_std:.4f}"

                # Compute speedup relative to scipy
                if ref_time and ref_time > 0:
                    speedup = ref_time / total_mean
                    speedup_str = f"{speedup:.2f}x"
                else:
                    speedup_str = "---"

                print(f"  {backend_name:<10} {total_str:<20} {eig_str:<25} {speedup_str:<10}")
            else:
                print(f"  {backend_name:<10} {'FAILED':<20} {result.get('error', 'Unknown error'):<25}")

    # Print interpretation
    print("\n" + "="*70)
    print("INTERPRETATION")
    print("="*70)
    print("""
For the tested configurations (iterations 3-5, fraction 0.7):
  - System sizes: 44 to 716 nodes
  - These are VERY SMALL systems

Expected behavior:
  - scipy/numpy will likely be FASTER than cupy at these sizes
  - cupy has GPU overhead that dominates for small problems
  - cupy typically becomes beneficial at ~5k-10k+ nodes

Recommendations:
  - For cluster deployment with large systems (iter 8+), use cupy
  - For quick tests with small systems, scipy is fine
  - Run with --extended flag to see cupy speedup at larger sizes

To test larger sizes, run:
  python benchmark_mcg_backends.py --extended
    """)


def main():
    """Run the benchmark."""
    import argparse

    parser = argparse.ArgumentParser(description="Benchmark MCG eigenvalue backends")
    parser.add_argument('--extended', action='store_true',
                        help='Include extended configurations (iter 6-8, larger sizes)')
    args = parser.parse_args()

    configs = TEST_CONFIGS.copy()

    if args.extended:
        print("\n⚠ Including EXTENDED configurations (this will take longer)")
        configs.extend(EXTENDED_CONFIGS)

    run_benchmark(configs, use_extended=args.extended)

    return 0


if __name__ == "__main__":
    sys.exit(main())
