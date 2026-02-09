#!/usr/bin/env python
"""
Quick test: Compare eigenvalue vs SLQ vs expm_multiply for entropy computation.

SLQ (Stochastic Lanczos Quadrature) is a middle-ground method:
- Avoids full eigendecomposition (good for large graphs)
- Uses Lanczos tridiagonalization (more efficient than expm_multiply)
- Should handle large tau ranges better
"""

import numpy as np
import time
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from lrgsglib.graphs.nx import MultiplicativeCascadeGraphNX as MultiplicativeCascadeGraph
from lrgsglib.utils.lrg.infocomm import (
    compute_entropy_observables_from_eigenvalues,
    compute_entropy_observables_slq,
    compute_entropy_observables_expm_multiply
)


def test_all_methods(iterations, fraction=0.7, p1=1.0, p2=0.8, p3=0.8, p4=0.9,
                     entropy_steps=50, t1=-2, t2=7):
    """Test all three entropy computation methods."""

    print(f"\n{'='*70}")
    print(f"Testing iterations={iterations}, fraction={fraction}")
    print(f"        p1={p1}, p2={p2}, p3={p3}, p4={p4}")
    print(f"        tau range: [10^{t1}, 10^{t2}], steps={entropy_steps}")
    print(f"{'='*70}")

    # Generate graph
    G = MultiplicativeCascadeGraph(
        p1=p1, p2=p2, p3=p3, p4=p4,
        fraction=fraction,
        iterations=iterations,
        pflip=0.1,
        seed=42
    )
    G.upd_graph_matrices()

    N = G.N
    E = G.gr['G'].number_of_edges()
    print(f"Graph: N={N:,} nodes, E={E:,} edges\n")

    results = {}

    # ========== Method 1: Eigenvalue ==========
    print("Method 1: Eigenvalue-based")
    try:
        t0 = time.time()
        G.compute_laplacian_spectrum(backend='scipy')
        t_spectrum = time.time() - t0

        t0 = time.time()
        S1, C1, V1, tau1 = compute_entropy_observables_from_eigenvalues(
            eigenvalues=G.eigv,
            num_nodes=G.N,
            steps=entropy_steps,
            t1=t1,
            t2=t2
        )
        t_entropy = time.time() - t0
        t_total = t_spectrum + t_entropy

        print(f"  Spectrum: {t_spectrum:.3f}s, Entropy: {t_entropy:.3f}s")
        print(f"  TOTAL: {t_total:.3f}s")

        results['eigenvalue'] = {
            'time': t_total,
            'S': S1,
            'C': C1,
            'tau': tau1,
            'success': True
        }
    except Exception as e:
        print(f"  FAILED: {e}")
        results['eigenvalue'] = {'success': False}

    # ========== Method 2: SLQ ==========
    print("\nMethod 2: SLQ (Stochastic Lanczos Quadrature)")
    try:
        t0 = time.time()
        S2, C2, V2, tau2 = compute_entropy_observables_slq(
            laplacian=G.slp,
            num_nodes=G.N,
            steps=entropy_steps,
            t1=t1,
            t2=t2,
            lanczos_steps=50,  # Lanczos iterations
            num_samples=100,   # Hutchinson samples
            seed=42
        )
        t_total = time.time() - t0

        print(f"  TOTAL: {t_total:.3f}s")

        results['slq'] = {
            'time': t_total,
            'S': S2,
            'C': C2,
            'tau': tau2,
            'success': True
        }

        # Compare with eigenvalue if available
        if results['eigenvalue']['success']:
            error = np.nanmean(np.abs(S1 - S2))
            max_error = np.nanmax(np.abs(S1 - S2))
            print(f"  Mean error: {error:.4f}, Max error: {max_error:.4f}")
            results['slq']['mean_error'] = error
            results['slq']['max_error'] = max_error

    except Exception as e:
        print(f"  FAILED: {e}")
        results['slq'] = {'success': False}

    # ========== Method 3: expm_multiply ==========
    print("\nMethod 3: expm_multiply")
    try:
        t0 = time.time()
        S3, C3, V3, tau3 = compute_entropy_observables_expm_multiply(
            L=G.slp,
            num_nodes=G.N,
            steps=entropy_steps,
            t1=t1,
            t2=t2,
            num_samples=100,
            seed=42,
            verbose=False
        )
        t_total = time.time() - t0

        print(f"  TOTAL: {t_total:.3f}s")

        results['expm'] = {
            'time': t_total,
            'S': S3,
            'C': C3,
            'tau': tau3,
            'success': True
        }

        # Compare with eigenvalue if available
        if results['eigenvalue']['success']:
            error = np.nanmean(np.abs(S1 - S3))
            max_error = np.nanmax(np.abs(S1 - S3))
            print(f"  Mean error: {error:.4f}, Max error: {max_error:.4f}")
            results['expm']['mean_error'] = error
            results['expm']['max_error'] = max_error

    except Exception as e:
        print(f"  FAILED: {e}")
        results['expm'] = {'success': False}

    # ========== Summary ==========
    print(f"\n{'='*70}")
    print("SUMMARY:")
    print(f"{'Method':<20} {'Time':>12} {'vs Eigenvalue':>15} {'Error':>12}")
    print("-"*70)

    eig_time = results['eigenvalue']['time'] if results['eigenvalue']['success'] else None

    for method_name, method_key in [('Eigenvalue', 'eigenvalue'),
                                     ('SLQ', 'slq'),
                                     ('expm_multiply', 'expm')]:
        if results[method_key]['success']:
            t = results[method_key]['time']
            if eig_time and method_key != 'eigenvalue':
                speedup = f"{eig_time/t:.2f}×"
                error = results[method_key].get('mean_error', 0)
                error_str = f"{error:.4f}"
            else:
                speedup = "baseline"
                error_str = "—"
            print(f"{method_name:<20} {t:>11.3f}s {speedup:>14} {error_str:>12}")
        else:
            print(f"{method_name:<20} {'FAILED':>12}")

    print(f"{'='*70}\n")

    return results


def main():
    print("\n" + "="*70)
    print("ENTROPY COMPUTATION: 3-WAY COMPARISON")
    print("Methods: Eigenvalue | SLQ | expm_multiply")
    print("="*70)

    # Test on progressively larger graphs
    test_sizes = [
        (3, "Small (N~150-200)"),
        (4, "Medium (N~600-700)"),
        (5, "Large (N~2500-3000)"),
    ]

    all_results = []

    for iterations, description in test_sizes:
        print(f"\n### {description} ###")
        result = test_all_methods(iterations, entropy_steps=50, t1=-2, t2=7)
        result['iterations'] = iterations
        all_results.append(result)

    # Final comparison table
    print("\n" + "="*70)
    print("FINAL COMPARISON")
    print("="*70)
    print(f"{'N':>8} {'Eigenvalue':>12} {'SLQ':>12} {'expm_multiply':>14} {'Winner':>12}")
    print("-"*70)

    for res in all_results:
        N = res['eigenvalue']['S'].shape[0] if res['eigenvalue']['success'] else "?"

        times = {}
        for key in ['eigenvalue', 'slq', 'expm']:
            if res[key]['success']:
                times[key] = res[key]['time']

        if times:
            winner = min(times.keys(), key=lambda k: times[k])
            winner_name = {'eigenvalue': 'Eigenvalue', 'slq': 'SLQ', 'expm': 'expm'}[winner]
        else:
            winner_name = "?"

        eig_str = f"{times.get('eigenvalue', 0):.2f}s" if 'eigenvalue' in times else "FAIL"
        slq_str = f"{times.get('slq', 0):.2f}s" if 'slq' in times else "FAIL"
        expm_str = f"{times.get('expm', 0):.2f}s" if 'expm' in times else "FAIL"

        print(f"{N:>8} {eig_str:>12} {slq_str:>12} {expm_str:>14} {winner_name:>12}")

    print("\nConclusion:")
    print("  - Eigenvalue: Best for N < 5000, requires full spectrum")
    print("  - SLQ: Middle ground, may be faster for large tau ranges")
    print("  - expm_multiply: Slowest, only use when eigendecomp impossible")


if __name__ == "__main__":
    main()
