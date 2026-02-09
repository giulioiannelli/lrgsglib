#!/usr/bin/env python
"""
Benchmark: Compare eigenvalue-based vs expm_multiply entropy computation.

Measures timing for:
1. Graph generation
2. Spectral computation (eigenvalues)
3. Entropy computation (eigenvalue-based vs expm_multiply)

Usage:
    # Quick benchmark (N < 1000)
    python bench_entropy_methods.py --quick

    # Full benchmark (N up to 10000+, can be slow!)
    python bench_entropy_methods.py --long

    # Custom parameters
    python bench_entropy_methods.py --iterations 4 --fraction 0.4
"""

import argparse
import time
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from lrgsglib.graphs.nx import MultiplicativeCascadeGraphNX as MultiplicativeCascadeGraph
from lrgsglib.utils.lrg.infocomm import (
    compute_entropy_observables_from_eigenvalues,
    compute_entropy_observables_expm_multiply
)

try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


def benchmark_single_graph(iterations, fraction, pflip, entropy_steps, entropy_range,
                          p1=1.0, p2=0.8, p3=0.8, p4=0.9, num_samples=150):
    """Benchmark both methods on a single graph."""
    print(f"\n{'='*70}")
    print(f"Benchmark: iterations={iterations}, fraction={fraction}, pflip={pflip}")
    print(f"           p1={p1}, p2={p2}, p3={p3}, p4={p4}")
    print(f"{'='*70}")

    # ========== 1. Graph Generation ==========
    t0 = time.time()
    G = MultiplicativeCascadeGraph(
        p1=p1, p2=p2, p3=p3, p4=p4,
        fraction=fraction,
        iterations=iterations,
        pflip=pflip,
        seed=42
    )
    t_gen = time.time() - t0

    N = G.N
    E = G.gr['G'].number_of_edges()
    print(f"Graph: N={N:,} nodes, E={E:,} edges")
    print(f"  Generation time: {t_gen:.3f}s")

    # Compute matrices
    G.upd_graph_matrices()

    results = {
        'N': N,
        'E': E,
        't_generation': t_gen
    }

    # Initialize data storage
    S1 = C1 = tau1 = None
    S2 = C2 = tau2 = None

    # ========== 2. Eigenvalue Method ==========
    print(f"\nMethod 1: Eigenvalue-based")
    try:
        t0 = time.time()
        G.compute_laplacian_spectrum(backend='scipy')
        t_spectrum = time.time() - t0
        print(f"  Spectrum computation: {t_spectrum:.3f}s")

        t0 = time.time()
        S1, C1, V1, tau1 = compute_entropy_observables_from_eigenvalues(
            eigenvalues=G.eigv,
            num_nodes=G.N,
            steps=entropy_steps,
            t1=entropy_range[0],
            t2=entropy_range[1]
        )
        t_entropy1 = time.time() - t0
        print(f"  Entropy computation: {t_entropy1:.3f}s")
        print(f"  TOTAL: {t_spectrum + t_entropy1:.3f}s")

        results['eigenvalue_spectrum_time'] = t_spectrum
        results['eigenvalue_entropy_time'] = t_entropy1
        results['eigenvalue_total_time'] = t_spectrum + t_entropy1
        results['eigenvalue_success'] = True

    except Exception as e:
        print(f"  FAILED: {e}")
        results['eigenvalue_success'] = False

    # ========== 3. expm_multiply Method ==========
    print(f"\nMethod 2: expm_multiply")
    try:
        t0 = time.time()
        S2, C2, V2, tau2 = compute_entropy_observables_expm_multiply(
            L=G.slp,  # Signed Laplacian (D_s - A)
            num_nodes=G.N,
            steps=entropy_steps,
            t1=entropy_range[0],
            t2=entropy_range[1],
            num_samples=num_samples,
            seed=42,
            verbose=False
        )
        t_expm = time.time() - t0
        print(f"  Entropy computation: {t_expm:.3f}s")
        print(f"  TOTAL: {t_expm:.3f}s (no spectrum needed!)")

        results['expm_time'] = t_expm
        results['expm_success'] = True

        # Compare results if both succeeded
        if results.get('eigenvalue_success'):
            error = np.nanmean(np.abs(S1 - S2))
            print(f"\n  Mean error: {error:.2e}")
            results['mean_error'] = error

    except Exception as e:
        print(f"  FAILED: {e}")
        results['expm_success'] = False

    # Store data for plotting
    results['S1'] = S1
    results['C1'] = C1
    results['tau1'] = tau1
    results['S2'] = S2
    results['C2'] = C2
    results['tau2'] = tau2

    # ========== 4. Summary ==========
    print(f"\n{'='*70}")
    print("Summary:")
    if results.get('eigenvalue_success') and results.get('expm_success'):
        speedup = results['eigenvalue_total_time'] / results['expm_time']
        print(f"  Eigenvalue method: {results['eigenvalue_total_time']:.3f}s")
        print(f"  expm_multiply method: {results['expm_time']:.3f}s")
        print(f"  Speedup: {speedup:.2f}×")
        results['speedup'] = speedup
    else:
        if not results.get('eigenvalue_success'):
            print("  Eigenvalue method: FAILED")
        if not results.get('expm_success'):
            print("  expm_multiply method: FAILED")

    return results


def plot_comparison(result, output_dir):
    """Generate comparison plots for entropy and specific heat."""
    if not MATPLOTLIB_AVAILABLE:
        print("  Matplotlib not available, skipping plots")
        return None

    if not result.get('eigenvalue_success') or not result.get('expm_success'):
        print("  One method failed, skipping plots")
        return None

    S1 = result['S1']
    C1 = result['C1']
    tau1 = result['tau1']
    S2 = result['S2']
    C2 = result['C2']
    tau2 = result['tau2']
    N = result['N']

    if S1 is None or S2 is None:
        return None

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Plot entropy
    ax1.plot(tau1, S1, 'b-', linewidth=2, label='Eigenvalue method', alpha=0.8)
    ax1.plot(tau2, S2, 'r--', linewidth=2, label='expm_multiply', alpha=0.6)
    ax1.set_xscale('log')
    ax1.set_xlabel(r'$\tau$', fontsize=12)
    ax1.set_ylabel(r'$S(\tau)$', fontsize=12)
    ax1.set_title(f'Entropy Comparison (N={N:,} nodes)', fontsize=14, fontweight='bold')
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)

    # Plot specific heat
    ax2.plot(tau1, C1, 'b-', linewidth=2, label='Eigenvalue method', alpha=0.8)
    ax2.plot(tau2, C2, 'r--', linewidth=2, label='expm_multiply', alpha=0.6)
    ax2.set_xscale('log')
    ax2.set_xlabel(r'$\tau$', fontsize=12)
    ax2.set_ylabel(r'$C(\tau)$', fontsize=12)
    ax2.set_title(f'Specific Heat Comparison (N={N:,} nodes)', fontsize=14, fontweight='bold')
    ax2.legend(loc='best')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    # Save plot
    plot_file = output_dir / f"entropy_comparison_N{N:06d}.png"
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f"  Plot saved: {plot_file}")
    return plot_file


def main():
    parser = argparse.ArgumentParser(description="Benchmark entropy computation methods")
    parser.add_argument('--quick', action='store_true',
                       help='Quick benchmark (small graphs, < 30s)')
    parser.add_argument('--long', action='store_true',
                       help='Long benchmark (large graphs N~1000+, can take minutes)')
    parser.add_argument('--iterations', type=int, nargs='+',
                       help='Graph iterations to test (overrides presets)')
    parser.add_argument('--fraction', type=float, default=None,
                       help='Node sampling fraction (overrides mode default)')
    parser.add_argument('--pflip', type=float, default=0.1,
                       help='Edge flip probability')
    parser.add_argument('--entropy-steps', type=int, default=50,
                       help='Number of entropy time points')
    parser.add_argument('--entropy-range', type=int, nargs=2, default=[-2, 3],
                       help='Entropy time range (t1 t2), e.g., [-2, 3] = [10^-2, 10^3]')
    parser.add_argument('--num-samples', type=int, default=150,
                       help='Number of Hutchinson samples for expm_multiply (default: 150)')
    parser.add_argument('--p1', type=float, default=1.0,
                       help='MCG probability p1 (default: 1.0)')
    parser.add_argument('--p2', type=float, default=0.8,
                       help='MCG probability p2 (default: 0.8)')
    parser.add_argument('--p3', type=float, default=0.8,
                       help='MCG probability p3 (default: 0.8)')
    parser.add_argument('--p4', type=float, default=0.9,
                       help='MCG probability p4 (default: 0.9)')
    args = parser.parse_args()

    # Set up output directory
    TEST_TMP_DIR = Path(__file__).parent / ".tmp"
    TEST_TMP_DIR.mkdir(exist_ok=True)

    # Determine test configurations
    if args.iterations:
        iterations_list = args.iterations
        fraction = args.fraction if args.fraction is not None else 0.7
        mode_desc = "custom"
    elif args.quick:
        iterations_list = [2, 3]  # With fraction=0.7: N ≈ 10-20 nodes (< 30s total)
        fraction = args.fraction if args.fraction is not None else 0.7
        mode_desc = "QUICK MODE: small graphs only"
    elif args.long:
        # With fraction=0.7: iterations 4,5,6 → N~50-300
        iterations_list = [4, 5, 6]
        fraction = args.fraction if args.fraction is not None else 0.7
        mode_desc = "LONG MODE: medium graphs (N~50-300), finds expm_multiply crossover"
    else:
        print("\nUsage: python bench_entropy_methods.py --quick   OR   --long")
        print("  --quick: Fast benchmark (< 30s, small graphs)")
        print("  --long:  Medium graphs (N~50-300, finds expm_multiply crossover)")
        return 1

    print("\n" + "="*70)
    print("ENTROPY METHODS PERFORMANCE BENCHMARK")
    print("="*70)
    print(f"Config: fraction={fraction}, pflip={args.pflip}")
    print(f"        p1={args.p1}, p2={args.p2}, p3={args.p3}, p4={args.p4}")
    print(f"        entropy_steps={args.entropy_steps}, range=[10^{args.entropy_range[0]}, 10^{args.entropy_range[1]}]")
    print(f"        num_samples={args.num_samples} (signed Laplacian)")
    print(f"\n[{mode_desc}]")

    # Run benchmarks
    all_results = []
    for it in iterations_list:
        result = benchmark_single_graph(
            iterations=it,
            fraction=fraction,
            pflip=args.pflip,
            entropy_steps=args.entropy_steps,
            entropy_range=args.entropy_range,
            p1=args.p1,
            p2=args.p2,
            p3=args.p3,
            p4=args.p4,
            num_samples=args.num_samples
        )
        all_results.append(result)

        # Generate plot for this result
        if MATPLOTLIB_AVAILABLE:
            plot_comparison(result, TEST_TMP_DIR)

    # Final summary
    print("\n" + "="*70)
    print("FINAL SUMMARY")
    print("="*70)
    print(f"{'N':>10} {'Eigenvalue':>12} {'expm_multiply':>14} {'Speedup':>10}")
    print("-"*70)
    for r in all_results:
        if r.get('speedup'):
            print(f"{r['N']:>10,} {r['eigenvalue_total_time']:>11.3f}s {r['expm_time']:>13.3f}s {r['speedup']:>9.2f}×")
        else:
            status_eig = "OK" if r.get('eigenvalue_success') else "FAIL"
            status_expm = "OK" if r.get('expm_success') else "FAIL"
            print(f"{r['N']:>10,} {status_eig:>12} {status_expm:>14} {'N/A':>10}")

    print("\nNote: Speedup > 1.0 means expm_multiply is faster")
    print("      Crossover typically at N ≈ 5,000-10,000 nodes")

    # Save results
    import pickle
    output_file = TEST_TMP_DIR / "bench_entropy_methods.pkl"
    with open(output_file, 'wb') as f:
        pickle.dump({
            'config': vars(args),
            'fraction': fraction,
            'results': all_results
        }, f)
    print(f"\nResults saved to: {output_file}")

    if MATPLOTLIB_AVAILABLE:
        print(f"Plots saved to: {TEST_TMP_DIR}/")
    else:
        print("(Matplotlib not available, no plots generated)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
