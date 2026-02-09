#!/usr/bin/env python
"""
Comprehensive Entropy Computation Benchmark

Tests ALL methods for computing S[ρ](τ):
1. Eigenvalue-based (exact):
   - Dense scipy (eigh)
   - Sparse scipy (eigsh, gets N-2 eigenvalues)
   - Dense cupy (GPU, if available)
2. SLQ (Stochastic Lanczos Quadrature) - various parameter regimes
3. expm_multiply (matrix exponential approximation)

Each method has 60-second timeout. Compares accuracy, speed, feasibility.
Visual comparison of exact vs approximate methods when both succeed.

Usage:
    python bench_entropy_comprehensive.py
"""

import numpy as np
import time
import signal
import warnings
from pathlib import Path
import sys
from typing import Optional, Tuple, Dict, List
from contextlib import contextmanager

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from lrgsglib.graphs.nx import MultiplicativeCascadeGraphNX as MultiplicativeCascadeGraph
from lrgsglib.utils.lrg.infocomm import (
    compute_entropy_observables_from_eigenvalues,
    compute_entropy_observables_slq,
    compute_entropy_observables_expm_multiply
)

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

# Global timeout flag
class TimeoutException(Exception):
    pass

def timeout_handler(signum, frame):
    raise TimeoutException("Computation exceeded 60 second timeout")

@contextmanager
def time_limit(seconds):
    """Context manager for timeout."""
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)


class EntropyBenchmark:
    """Comprehensive entropy computation benchmark."""

    def __init__(self, timeout_seconds=60):
        self.timeout = timeout_seconds
        self.results = []

    def test_eigenvalue_dense_scipy(self, G, steps, t1, t2):
        """Test dense scipy eigendecomposition (full spectrum)."""
        method_name = "Eigenvalue (dense scipy)"
        try:
            with time_limit(self.timeout):
                t0 = time.time()
                G.compute_laplacian_spectrum(backend='scipy')
                t_spec = time.time() - t0

                S, C, V, tau = compute_entropy_observables_from_eigenvalues(
                    eigenvalues=G.eigv,
                    num_nodes=G.N,
                    steps=steps,
                    t1=t1,
                    t2=t2
                )
                t_total = time.time() - t0

                return {
                    'method': method_name,
                    'success': True,
                    'time': t_total,
                    'time_spectrum': t_spec,
                    'S': S,
                    'C': C,
                    'tau': tau,
                    'error': None,
                    'status': 'OK'
                }
        except TimeoutException:
            return {
                'method': method_name,
                'success': False,
                'time': self.timeout,
                'status': f'TIMEOUT (>{self.timeout}s)',
                'error': 'Computation too slow'
            }
        except MemoryError:
            return {
                'method': method_name,
                'success': False,
                'status': 'MEMORY ERROR',
                'error': f'Dense matrix too large for N={G.N}'
            }
        except Exception as e:
            return {
                'method': method_name,
                'success': False,
                'status': 'ERROR',
                'error': str(e)
            }

    def test_eigenvalue_sparse_scipy(self, G, steps, t1, t2):
        """Test sparse scipy eigendecomposition (N-2 eigenvalues)."""
        method_name = "Eigenvalue (sparse scipy)"
        try:
            with time_limit(self.timeout):
                t0 = time.time()
                # Force sparse mode - gets N-2 eigenvalues
                G.compute_laplacian_spectrum(backend='scipy', keep_sparse=True)
                t_spec = time.time() - t0

                if G.eigv is None or len(G.eigv) < G.N - 10:
                    raise ValueError(f"Sparse mode failed: got {len(G.eigv) if G.eigv is not None else 0} eigenvalues, expected ~{G.N-2}")

                S, C, V, tau = compute_entropy_observables_from_eigenvalues(
                    eigenvalues=G.eigv,
                    num_nodes=G.N,
                    steps=steps,
                    t1=t1,
                    t2=t2
                )
                t_total = time.time() - t0

                return {
                    'method': method_name,
                    'success': True,
                    'time': t_total,
                    'time_spectrum': t_spec,
                    'S': S,
                    'C': C,
                    'tau': tau,
                    'error': None,
                    'status': f'OK (N-2={len(G.eigv)} eigenvalues)',
                    'num_eigenvalues': len(G.eigv)
                }
        except TimeoutException:
            return {
                'method': method_name,
                'success': False,
                'time': self.timeout,
                'status': f'TIMEOUT (>{self.timeout}s)',
                'error': 'Computation too slow'
            }
        except Exception as e:
            return {
                'method': method_name,
                'success': False,
                'status': 'ERROR',
                'error': str(e)
            }

    def test_eigenvalue_gpu_cupy(self, G, steps, t1, t2):
        """Test GPU eigendecomposition with cupy."""
        method_name = "Eigenvalue (dense cupy GPU)"
        try:
            import cupy
            with time_limit(self.timeout):
                t0 = time.time()
                G.compute_laplacian_spectrum(backend='cupy')
                t_spec = time.time() - t0

                S, C, V, tau = compute_entropy_observables_from_eigenvalues(
                    eigenvalues=G.eigv,
                    num_nodes=G.N,
                    steps=steps,
                    t1=t1,
                    t2=t2
                )
                t_total = time.time() - t0

                return {
                    'method': method_name,
                    'success': True,
                    'time': t_total,
                    'time_spectrum': t_spec,
                    'S': S,
                    'C': C,
                    'tau': tau,
                    'error': None,
                    'status': 'OK'
                }
        except ImportError:
            return {
                'method': method_name,
                'success': False,
                'status': 'UNAVAILABLE',
                'error': 'cupy not installed'
            }
        except TimeoutException:
            return {
                'method': method_name,
                'success': False,
                'time': self.timeout,
                'status': f'TIMEOUT (>{self.timeout}s)',
                'error': 'Computation too slow'
            }
        except Exception as e:
            # Common: CUDA out of memory
            return {
                'method': method_name,
                'success': False,
                'status': 'ERROR',
                'error': str(e)[:100]  # Truncate long CUDA errors
            }

    def test_slq(self, G, steps, t1, t2, lanczos_steps, num_samples):
        """Test SLQ with specific parameters."""
        method_name = f"SLQ (L={lanczos_steps}, S={num_samples})"
        try:
            with time_limit(self.timeout):
                t0 = time.time()
                S, C, V, tau = compute_entropy_observables_slq(
                    laplacian=G.slp,
                    num_nodes=G.N,
                    steps=steps,
                    t1=t1,
                    t2=t2,
                    lanczos_steps=lanczos_steps,
                    num_samples=num_samples,
                    seed=42
                )
                t_total = time.time() - t0

                # Check for numerical issues (entropy > 1 means unstable)
                max_entropy = np.nanmax(S)
                if max_entropy > 1.5:
                    warning = f"Unstable (max S={max_entropy:.2f} > 1)"
                else:
                    warning = None

                return {
                    'method': method_name,
                    'success': True,
                    'time': t_total,
                    'S': S,
                    'C': C,
                    'tau': tau,
                    'error': None,
                    'status': 'OK' if warning is None else warning,
                    'lanczos_steps': lanczos_steps,
                    'num_samples': num_samples
                }
        except TimeoutException:
            return {
                'method': method_name,
                'success': False,
                'time': self.timeout,
                'status': f'TIMEOUT (>{self.timeout}s)',
                'error': 'Computation too slow'
            }
        except Exception as e:
            return {
                'method': method_name,
                'success': False,
                'status': 'ERROR',
                'error': str(e)
            }

    def test_expm_multiply(self, G, steps, t1, t2, num_samples):
        """Test expm_multiply method."""
        method_name = f"expm_multiply (S={num_samples})"
        try:
            with time_limit(self.timeout):
                t0 = time.time()
                S, C, V, tau = compute_entropy_observables_expm_multiply(
                    L=G.slp,
                    num_nodes=G.N,
                    steps=steps,
                    t1=t1,
                    t2=t2,
                    num_samples=num_samples,
                    seed=42,
                    verbose=False
                )
                t_total = time.time() - t0

                return {
                    'method': method_name,
                    'success': True,
                    'time': t_total,
                    'S': S,
                    'C': C,
                    'tau': tau,
                    'error': None,
                    'status': 'OK',
                    'num_samples': num_samples
                }
        except TimeoutException:
            return {
                'method': method_name,
                'success': False,
                'time': self.timeout,
                'status': f'TIMEOUT (>{self.timeout}s)',
                'error': 'Computation too slow (large tau)'
            }
        except Exception as e:
            return {
                'method': method_name,
                'success': False,
                'status': 'ERROR',
                'error': str(e)
            }

    def benchmark_graph(self, p1, p2, p3, p4, iterations, fraction, steps=50, t1=-2, t2=7,
                        skip_eigenvalue=False, skip_expm=False, only_slq=False):
        """Run all methods on a single graph configuration."""
        print(f"\n{'='*80}")
        print(f"Graph: p1={p1}, p2={p2}, p3={p3}, p4={p4}, iter={iterations}, frac={fraction}")
        print(f"Tau range: [10^{t1}, 10^{t2}], steps={steps}")
        print(f"{'='*80}")

        # Generate graph
        t0 = time.time()
        G = MultiplicativeCascadeGraph(
            p1=p1, p2=p2, p3=p3, p4=p4,
            fraction=fraction,
            iterations=iterations,
            pflip=0.1,
            seed=42
        )
        G.upd_graph_matrices()
        t_gen = time.time() - t0

        N = G.N
        E = G.gr['G'].number_of_edges()
        print(f"Generated: N={N:,} nodes, E={E:,} edges ({t_gen:.2f}s)")
        print()

        results = {
            'N': N,
            'E': E,
            'p1': p1, 'p2': p2, 'p3': p3, 'p4': p4,
            'iterations': iterations,
            'fraction': fraction,
            't_gen': t_gen,
            'methods': []
        }

        # Test all eigenvalue methods
        if not skip_eigenvalue and not only_slq:
            print("Testing eigenvalue-based methods...")
            eigenvalue_tests = [
                (self.test_eigenvalue_dense_scipy, "Dense scipy"),
                (self.test_eigenvalue_sparse_scipy, "Sparse scipy"),
                (self.test_eigenvalue_gpu_cupy, "GPU cupy"),
            ]

            for test_func, desc in eigenvalue_tests:
                # Skip GPU for large graphs to avoid timeout -> fallback hang
                if "GPU" in desc and N > 5000:
                    print(f"  {desc}... ⊘ SKIPPED (N > 5000)")
                    results['methods'].append({
                        'method': desc,
                        'success': False,
                        'status': 'SKIPPED (N > 5000)',
                        'error': 'Graph too large for GPU timeout handling'
                    })
                    continue

                print(f"  {desc}...", end=' ', flush=True)
                result = test_func(G, steps, t1, t2)
                results['methods'].append(result)
                if result['success']:
                    print(f"✓ {result['time']:.2f}s")
                else:
                    print(f"✗ {result['status']}")
        else:
            print("Eigenvalue methods... ⊘ SKIPPED (--skip-eigenvalue or --only-slq)")

        # Find exact reference (prefer dense scipy, fallback to sparse)
        exact_result = None
        for r in results['methods']:
            if r['success'] and 'eigenvalue' in r['method'].lower():
                exact_result = r
                break

        # Test SLQ with various parameter regimes
        print("\nTesting SLQ methods...")
        slq_configs = [
            (50, 50, "Fast/Low accuracy"),
            (100, 100, "Balanced"),
            (200, 100, "High accuracy/Medium"),
            (300, 100, "High accuracy/Slow"),
        ]

        for lanczos, samples, desc in slq_configs:
            print(f"  SLQ (L={lanczos}, S={samples}) [{desc}]...", end=' ', flush=True)
            result = self.test_slq(G, steps, t1, t2, lanczos, samples)

            # Compute error vs exact if available
            if result['success'] and exact_result is not None:
                error = np.nanmean(np.abs(exact_result['S'] - result['S']))
                max_error = np.nanmax(np.abs(exact_result['S'] - result['S']))
                result['mean_error'] = error
                result['max_error'] = max_error
                print(f"✓ {result['time']:.2f}s (err={error:.4f})")
            elif result['success']:
                print(f"✓ {result['time']:.2f}s (no reference)")
            else:
                print(f"✗ {result['status']}")

            results['methods'].append(result)

        # Test expm_multiply
        if not skip_expm and not only_slq:
            print("\nTesting expm_multiply method...")
            print(f"  expm_multiply (S=100)...", end=' ', flush=True)
            result = self.test_expm_multiply(G, steps, t1, t2, num_samples=100)

            if result['success'] and exact_result is not None:
                error = np.nanmean(np.abs(exact_result['S'] - result['S']))
                max_error = np.nanmax(np.abs(exact_result['S'] - result['S']))
                result['mean_error'] = error
                result['max_error'] = max_error
                print(f"✓ {result['time']:.2f}s (err={error:.4f})")
            elif result['success']:
                print(f"✓ {result['time']:.2f}s (no reference)")
            else:
                print(f"✗ {result['status']}")

            results['methods'].append(result)
        else:
            print("\nexpm_multiply method... ⊘ SKIPPED (--skip-expm or --only-slq)")

        # Summary table
        print(f"\n{'='*80}")
        print("SUMMARY")
        print(f"{'='*80}")
        print(f"{'Method':<35} {'Time':>12} {'Status':>15} {'Error':>12}")
        print("-"*80)

        for r in results['methods']:
            method = r['method']
            if r['success']:
                time_str = f"{r['time']:.2f}s"
                status = r.get('status', 'OK')
                error_str = f"{r.get('mean_error', 0):.4f}" if 'mean_error' in r else "—"
            else:
                time_str = "—"
                status = r.get('status', 'FAILED')
                error_str = "—"

            print(f"{method:<35} {time_str:>12} {status:>15} {error_str:>12}")

        # Generate comparison plot if we have exact + approximations
        if MATPLOTLIB_AVAILABLE and exact_result is not None:
            self.plot_comparison(results, exact_result)
        elif MATPLOTLIB_AVAILABLE:
            # Plot SLQ-only comparison if no exact reference
            slq_results = [r for r in results['methods'] if 'SLQ' in r['method'] and r['success']]
            if len(slq_results) >= 2:
                self.plot_slq_comparison(results, slq_results)

        return results

    def plot_comparison(self, results, exact_result):
        """Generate visual comparison plot."""
        N = results['N']
        output_dir = Path(__file__).parent / ".tmp"
        output_dir.mkdir(exist_ok=True)

        # Collect successful approximate methods
        approx_methods = [r for r in results['methods']
                         if r['success'] and r != exact_result and 'mean_error' in r]

        if not approx_methods:
            return

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))

        # Plot exact
        tau = exact_result['tau']
        ax1.plot(tau, exact_result['S'], 'k-', linewidth=3, label='Exact (eigenvalues)', alpha=0.8, zorder=10)
        ax2.plot(tau, exact_result['C'], 'k-', linewidth=3, label='Exact (eigenvalues)', alpha=0.8, zorder=10)

        # Plot approximations
        colors = plt.cm.tab10(np.linspace(0, 1, len(approx_methods)))
        for i, r in enumerate(approx_methods):
            label = f"{r['method']} (err={r.get('mean_error', 0):.3f})"
            ax1.plot(r['tau'], r['S'], '--', linewidth=2, label=label, alpha=0.6, color=colors[i])
            ax2.plot(r['tau'], r['C'], '--', linewidth=2, label=label, alpha=0.6, color=colors[i])

        # Formatting
        ax1.set_xscale('log')
        ax1.set_xlabel(r'$\tau$', fontsize=12)
        ax1.set_ylabel(r'$S(\tau)$', fontsize=12)
        ax1.set_title(f'Entropy Comparison (N={N:,} nodes)', fontsize=14, fontweight='bold')
        ax1.legend(loc='best', fontsize=9)
        ax1.grid(True, alpha=0.3)

        ax2.set_xscale('log')
        ax2.set_xlabel(r'$\tau$', fontsize=12)
        ax2.set_ylabel(r'$C(\tau)$', fontsize=12)
        ax2.set_title(f'Specific Heat Comparison (N={N:,} nodes)', fontsize=14, fontweight='bold')
        ax2.legend(loc='best', fontsize=9)
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        # Save
        plot_file = output_dir / f"entropy_comparison_comprehensive_N{N:06d}.png"
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close(fig)

        print(f"\n  Plot saved: {plot_file}")

    def plot_slq_comparison(self, results, slq_results):
        """Generate SLQ-only comparison plot (no exact reference needed)."""
        N = results['N']
        output_dir = Path(__file__).parent / ".tmp"
        output_dir.mkdir(exist_ok=True)

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 14))

        # Plot all SLQ configurations
        colors = plt.cm.viridis(np.linspace(0, 1, len(slq_results)))
        max_entropy_overall = 0.0

        for i, r in enumerate(slq_results):
            tau = r['tau']
            S = r['S']
            C = r['C']

            max_entropy = np.nanmax(S)
            max_entropy_overall = max(max_entropy_overall, max_entropy)

            label = f"{r['method']} (max S={max_entropy:.3f})"
            ax1.plot(tau, S, linewidth=2, label=label, alpha=0.8, color=colors[i])
            ax2.plot(tau, C, linewidth=2, label=label, alpha=0.8, color=colors[i])

            # Plot error indicator in third panel (difference from best config)
            if i > 0:
                diff = np.abs(S - slq_results[0]['S'])
                ax3.plot(tau, diff, linewidth=2, label=f"{r['method']} vs {slq_results[0]['method']}",
                        alpha=0.7, color=colors[i])

        # Formatting panel 1: Entropy
        ax1.set_xscale('log')
        ax1.set_xlabel(r'$\tau$', fontsize=12)
        ax1.set_ylabel(r'$S(\tau)$', fontsize=12)
        ax1.set_title(f'SLQ Entropy Comparison (N={N:,} nodes)', fontsize=14, fontweight='bold')
        ax1.axhline(y=1.0, color='red', linestyle='--', linewidth=1, alpha=0.5, label='S=1 limit')
        ax1.legend(loc='best', fontsize=9)
        ax1.grid(True, alpha=0.3)

        # Add warning if entropy exceeds 1
        if max_entropy_overall > 1.05:
            ax1.text(0.5, 0.95, f'⚠ WARNING: Max entropy = {max_entropy_overall:.3f} > 1.0 (unstable!)',
                    transform=ax1.transAxes, ha='center', va='top',
                    bbox=dict(boxstyle='round', facecolor='red', alpha=0.3),
                    fontsize=11, fontweight='bold')

        # Formatting panel 2: Specific heat
        ax2.set_xscale('log')
        ax2.set_xlabel(r'$\tau$', fontsize=12)
        ax2.set_ylabel(r'$C(\tau)$', fontsize=12)
        ax2.set_title(f'Specific Heat Comparison (N={N:,} nodes)', fontsize=14, fontweight='bold')
        ax2.legend(loc='best', fontsize=9)
        ax2.grid(True, alpha=0.3)

        # Formatting panel 3: Differences
        if len(slq_results) > 1:
            ax3.set_xscale('log')
            ax3.set_yscale('log')
            ax3.set_xlabel(r'$\tau$', fontsize=12)
            ax3.set_ylabel(r'$|S - S_{reference}|$', fontsize=12)
            ax3.set_title(f'Entropy Differences (Reference: {slq_results[0]["method"]})',
                         fontsize=14, fontweight='bold')
            ax3.legend(loc='best', fontsize=9)
            ax3.grid(True, alpha=0.3)
        else:
            fig.delaxes(ax3)

        plt.tight_layout()

        # Save
        plot_file = output_dir / f"entropy_slq_comparison_N{N:06d}.png"
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close(fig)

        print(f"\n  SLQ comparison plot saved: {plot_file}")


def main():
    """Run comprehensive benchmark suite."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Comprehensive entropy computation benchmark",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all tests
  python bench_entropy_comprehensive.py

  # Run only large graph tests (skip medium)
  python bench_entropy_comprehensive.py --skip-size 0

  # Run only SLQ methods (skip eigenvalue and expm_multiply)
  python bench_entropy_comprehensive.py --only-slq

  # Run with longer timeout (2 minutes)
  python bench_entropy_comprehensive.py --timeout 120

  # Skip slow expm_multiply method
  python bench_entropy_comprehensive.py --skip-expm

  # Run only large graphs from config 0, SLQ only
  python bench_entropy_comprehensive.py --config 0 --size 1 --only-slq
        """
    )
    parser.add_argument(
        '--config', type=int, default=None,
        help='Run only specific config (0=high, 1=medium, 2=very high connectivity)'
    )
    parser.add_argument(
        '--size', type=int, default=None,
        help='Run only specific size (0=medium/~2-3k, 1=large/~10-15k, 2=very large/~30-50k)'
    )
    parser.add_argument(
        '--skip-size', type=int, nargs='+', default=[],
        help='Skip specific sizes (e.g., --skip-size 0 to skip medium graphs)'
    )
    parser.add_argument(
        '--skip-config', type=int, nargs='+', default=[],
        help='Skip specific configs'
    )
    parser.add_argument(
        '--timeout', type=int, default=60,
        help='Timeout per method in seconds (default: 60)'
    )
    parser.add_argument(
        '--skip-expm', action='store_true',
        help='Skip expm_multiply method (very slow for large graphs)'
    )
    parser.add_argument(
        '--skip-eigenvalue', action='store_true',
        help='Skip all eigenvalue methods (dense/sparse/GPU)'
    )
    parser.add_argument(
        '--only-slq', action='store_true',
        help='Run ONLY SLQ methods (skip everything else)'
    )

    args = parser.parse_args()

    print("="*80)
    print("COMPREHENSIVE ENTROPY COMPUTATION BENCHMARK")
    print("="*80)
    print(f"Timeout: {args.timeout} seconds per method")

    # Show what we're testing
    methods = []
    if not args.skip_eigenvalue and not args.only_slq:
        methods.append("Eigenvalue (dense/sparse/GPU)")
    methods.append("SLQ (4 configs)")
    if not args.skip_expm and not args.only_slq:
        methods.append("expm_multiply")
    print(f"Tests: {', '.join(methods)}")
    print("="*80)

    benchmark = EntropyBenchmark(timeout_seconds=args.timeout)

    # Test matrix: different probability configurations
    test_configs = [
        # (p1, p2, p3, p4, description)
        (1.0, 0.7, 0.7, 0.9, "High connectivity"),
        (1.0, 0.5, 0.5, 0.9, "Medium connectivity"),
        (1.0, 0.8, 0.8, 0.9, "Very high connectivity"),
    ]

    # Graph sizes
    test_sizes = [
        (5, 0.7, "Medium (N~2-3k)"),
        (6, 0.7, "Large (N~10-15k)"),
        (7, 0.6, "Very large (N~30-50k)"),
    ]

    all_results = []

    for config_idx, (p1, p2, p3, p4, config_desc) in enumerate(test_configs):
        # Skip if requested
        if args.config is not None and config_idx != args.config:
            continue
        if config_idx in args.skip_config:
            print(f"\n{'#'*80}")
            print(f"# SKIPPED Configuration {config_idx}: {config_desc}")
            print(f"{'#'*80}")
            continue

        print(f"\n{'#'*80}")
        print(f"# Configuration: {config_desc}")
        print(f"# p1={p1}, p2={p2}, p3={p3}, p4={p4}")
        print(f"{'#'*80}")

        for size_idx, (iterations, fraction, size_desc) in enumerate(test_sizes):
            # Skip if requested
            if args.size is not None and size_idx != args.size:
                continue
            if size_idx in args.skip_size:
                print(f"\n## SKIPPED {size_desc} (iterations={iterations}, fraction={fraction})")
                continue
            print(f"\n## {size_desc} (iterations={iterations}, fraction={fraction})")

            result = benchmark.benchmark_graph(
                p1=p1, p2=p2, p3=p3, p4=p4,
                iterations=iterations,
                fraction=fraction,
                steps=50,
                t1=-2,
                t2=7,
                skip_eigenvalue=args.skip_eigenvalue,
                skip_expm=args.skip_expm,
                only_slq=args.only_slq
            )
            all_results.append(result)

    # Final summary across all configurations
    if all_results:
        print(f"\n{'='*80}")
        print("FINAL SUMMARY - METHOD FEASIBILITY")
        print(f"{'='*80}")
        print(f"{'Graph Size':<15} {'Dense Scipy':>15} {'Sparse Scipy':>15} {'GPU Cupy':>15} {'SLQ Best':>15} {'expm':>10}")
        print("-"*95)
    else:
        print("\nNo tests were run (all skipped).")
        return

    for res in all_results:
        N = res['N']
        methods = {r['method']: r for r in res['methods']}

        # Extract status for each method type
        dense_status = "✓" if any('dense scipy' in m.lower() and methods[m]['success'] for m in methods) else "✗"
        sparse_status = "✓" if any('sparse scipy' in m.lower() and methods[m]['success'] for m in methods) else "✗"
        gpu_status = "✓" if any('cupy' in m.lower() and methods[m]['success'] for m in methods) else "✗"

        # Best SLQ (lowest error if any succeeded)
        slq_methods = [r for r in res['methods'] if 'SLQ' in r['method'] and r['success']]
        if slq_methods:
            best_slq = min(slq_methods, key=lambda x: x.get('mean_error', float('inf')))
            slq_status = f"✓ (L={best_slq.get('lanczos_steps', '?')})"
        else:
            slq_status = "✗"

        expm_status = "✓" if any('expm' in m.lower() and methods[m]['success'] for m in methods) else "✗"

        print(f"N={N:<12,} {dense_status:>15} {sparse_status:>15} {gpu_status:>15} {slq_status:>15} {expm_status:>10}")

    print("\n" + "="*80)
    print("RECOMMENDATIONS")
    print("="*80)
    print("N < 1,000:   Dense scipy (fastest, exact)")
    print("N ~ 1-5k:    Sparse scipy OR SLQ (L=200-300)")
    print("N ~ 5-20k:   SLQ (L=200-300) OR sparse scipy")
    print("N > 20k:     SLQ only (eigendecomp infeasible)")
    print("\nFor large tau (>10^6): All methods have numerical issues")
    print("                        Use tau cutoff and limiting values")
    print("="*80)


if __name__ == "__main__":
    main()
