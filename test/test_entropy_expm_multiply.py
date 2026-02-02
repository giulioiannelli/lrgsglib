#!/usr/bin/env python
"""
Test entropy computation with expm_multiply approach.

Validates the new entropy mode implementation by comparing against
eigenvalue-based computation on small graphs.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add src to path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph
from lrgsglib.utils.lrg.infocomm import (
    compute_entropy_observables_from_eigenvalues,
    compute_entropy_observables_expm_multiply
)


def test_expm_multiply_vs_eigenvalue_small_graph():
    """
    Test 1: Compare expm_multiply vs eigenvalue method on small graph.

    Expected: Errors decrease as num_samples increases.
    """
    print("\n" + "="*70)
    print("Test 1: expm_multiply vs Eigenvalue Method (Small Graph)")
    print("="*70)

    # Create small test graph
    print("\nCreating test graph...")
    G = MultiplicativeCascadeGraph(
        p1=0.8, p2=0.6, p3=0.6, p4=0.8,
        fraction=0.3,
        iterations=3,
        pflip=0.1,
        seed=42
    )
    print(f"Graph: N={G.N} nodes, E={G.gr['G'].number_of_edges()} edges")

    # Method 1: Eigenvalue-based (ground truth)
    print("\nMethod 1: Computing via eigenvalues...")
    G.compute_laplacian_spectrum(backend='scipy')
    S1, C1, V1, tau1 = compute_entropy_observables_from_eigenvalues(
        eigenvalues=G.eigv,
        num_nodes=G.N,
        steps=100,
        t1=-2,
        t2=5
    )
    print(f"  Computed {len(S1)} time points")

    # Method 2: expm_multiply with different num_samples
    print("\nMethod 2: Computing via expm_multiply...")
    G.upd_graph_matrices()  # Ensure slp is available

    errors_by_samples = []
    samples_list = [10, 30, 100]

    for num_samples in samples_list:
        S2, C2, V2, tau2 = compute_entropy_observables_expm_multiply(
            L=G.slp,
            num_nodes=G.N,
            steps=100,
            t1=-2,
            t2=5,
            num_samples=num_samples,
            seed=42
        )

        # Compute errors
        entropy_error = np.abs(S1 - S2)
        mean_error = np.mean(entropy_error)
        max_error = np.max(entropy_error)

        errors_by_samples.append({
            'num_samples': num_samples,
            'mean_error': mean_error,
            'max_error': max_error,
            'S': S2,
            'C': C2
        })

        print(f"  num_samples={num_samples:3d}: mean_error={mean_error:.6f}, max_error={max_error:.6f}")

    # Verify error is reasonably small with sufficient samples
    mean_errors = [e['mean_error'] for e in errors_by_samples]
    final_error = mean_errors[-1]
    # With 100 samples, error should be below 0.05 for this small graph
    error_threshold = 0.05
    if final_error < error_threshold:
        print(f"\n✓ PASS: Final error ({final_error:.6f}) is below threshold ({error_threshold})")
    else:
        print(f"\n✗ WARNING: Final error ({final_error:.6f}) exceeds threshold ({error_threshold})")

    # Create visualization
    TEST_TMP_DIR = Path(__file__).parent / ".tmp"
    TEST_TMP_DIR.mkdir(exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: Entropy comparison
    ax = axes[0, 0]
    ax.semilogx(tau1, S1, 'k-', linewidth=2, label='Eigenvalue (ground truth)')
    for result in errors_by_samples:
        ax.semilogx(tau1, result['S'], 'o-', alpha=0.6,
                   label=f"expm_multiply (m={result['num_samples']})")
    ax.set_xlabel('τ')
    ax.set_ylabel('Normalized Entropy S(τ)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Entropy Profiles')

    # Plot 2: Specific heat comparison
    ax = axes[0, 1]
    ax.semilogx(tau1, C1, 'k-', linewidth=2, label='Eigenvalue (ground truth)')
    for result in errors_by_samples:
        ax.semilogx(tau1, result['C'], 'o-', alpha=0.6,
                   label=f"expm_multiply (m={result['num_samples']})")
    ax.set_xlabel('τ')
    ax.set_ylabel('Specific Heat C(τ)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Specific Heat Profiles')

    # Plot 3: Absolute error vs τ
    ax = axes[1, 0]
    for result in errors_by_samples:
        error = np.abs(S1 - result['S'])
        ax.semilogx(tau1, error, 'o-', alpha=0.6,
                   label=f"m={result['num_samples']}")
    ax.set_xlabel('τ')
    ax.set_ylabel('|S_exact - S_expm|')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Absolute Error in Entropy')

    # Plot 4: Error convergence
    ax = axes[1, 1]
    samples_vals = [e['num_samples'] for e in errors_by_samples]
    mean_errs = [e['mean_error'] for e in errors_by_samples]
    max_errs = [e['max_error'] for e in errors_by_samples]

    ax.loglog(samples_vals, mean_errs, 'o-', label='Mean error')
    ax.loglog(samples_vals, max_errs, 's-', label='Max error')

    # Expected O(1/√m) convergence
    m_theory = np.array(samples_vals)
    theory_line = mean_errs[0] * (samples_vals[0] / m_theory)**0.5
    ax.loglog(m_theory, theory_line, 'k--', alpha=0.5, label='O(1/√m) theory')

    ax.set_xlabel('Number of samples (m)')
    ax.set_ylabel('Error')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_title('Error Convergence')

    plt.tight_layout()
    output_path = TEST_TMP_DIR / "entropy_expm_multiply_validation.png"
    plt.savefig(output_path, dpi=150)
    print(f"\n✓ Saved validation plot to: {output_path}")
    plt.close()

    # Assert error is reasonably small with sufficient samples
    assert final_error < error_threshold, f"Error ({final_error:.6f}) should be below {error_threshold}"


def test_entropy_mode_basic():
    """
    Test 2: Basic functionality test for entropy observables.

    Ensures that the function runs without errors and produces sensible output.
    """
    print("\n" + "="*70)
    print("Test 2: Basic Functionality Test")
    print("="*70)

    # Create graph
    print("\nCreating test graph...")
    G = MultiplicativeCascadeGraph(
        p1=0.9, p2=0.5, p3=0.5, p4=0.9,
        fraction=0.4,
        iterations=4,
        pflip=0.0,
        seed=123
    )
    G.upd_graph_matrices()
    print(f"Graph: N={G.N} nodes, E={G.gr['G'].number_of_edges()} edges")

    # Compute entropy
    print("\nComputing entropy observables...")
    S, C, V, tau = compute_entropy_observables_expm_multiply(
        L=G.slp,
        num_nodes=G.N,
        steps=200,
        t1=-2,
        t2=5,
        num_samples=50,
        seed=42
    )

    # Validate outputs
    checks = []

    # Check shapes
    checks.append(("Shape consistency",
                  len(S) == len(C) == len(V) == len(tau) == 200))

    # Check no NaNs
    checks.append(("No NaNs in entropy", not np.any(np.isnan(S))))
    checks.append(("No NaNs in specific heat", not np.any(np.isnan(C))))

    # Check finite values
    checks.append(("All finite", np.all(np.isfinite(S)) and np.all(np.isfinite(C))))

    # Check entropy bounds (allow slight overshoot for numerical reasons)
    checks.append(("Entropy in [0, 1.5]", np.all((S >= 0) & (S <= 1.5))))

    # Check tau grid is log-spaced
    log_tau = np.log10(tau)
    tau_spacing = np.diff(log_tau)
    checks.append(("Tau log-spaced",
                  np.allclose(tau_spacing, tau_spacing[0], rtol=1e-10)))

    # Print results
    print("\nValidation checks:")
    all_passed = True
    for check_name, passed in checks:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"  {status}: {check_name}")
        all_passed = all_passed and passed

    if all_passed:
        print("\n✓ All basic functionality tests passed!")
    else:
        print("\n✗ Some tests failed!")

    assert all_passed, "Some basic functionality tests failed"


def test_deterministic_seeding():
    """
    Test 3: Verify that seeding produces reproducible results.
    """
    print("\n" + "="*70)
    print("Test 3: Deterministic Seeding Test")
    print("="*70)

    # Create graph
    G = MultiplicativeCascadeGraph(
        p1=0.8, p2=0.6, p3=0.6, p4=0.8,
        fraction=0.3,
        iterations=3,
        pflip=0.15,
        seed=456
    )
    G.upd_graph_matrices()
    print(f"Graph: N={G.N} nodes")

    # Compute twice with same seed
    print("\nComputing entropy with seed=100 (run 1)...")
    S1, C1, V1, tau1 = compute_entropy_observables_expm_multiply(
        L=G.slp, num_nodes=G.N, steps=100, num_samples=30, seed=100
    )

    print("Computing entropy with seed=100 (run 2)...")
    S2, C2, V2, tau2 = compute_entropy_observables_expm_multiply(
        L=G.slp, num_nodes=G.N, steps=100, num_samples=30, seed=100
    )

    # Check match (allow small tolerance for numerical noise)
    entropy_match = np.allclose(S1, S2, rtol=1e-10, atol=1e-10)
    heat_match = np.allclose(C1, C2, rtol=1e-10, atol=1e-10)

    print(f"\nResults:")
    print(f"  Entropy exact match: {entropy_match}")
    print(f"  Specific heat exact match: {heat_match}")

    if entropy_match and heat_match:
        print("\n✓ PASS: Seeding produces deterministic results")
    else:
        print("\n✗ FAIL: Seeding did not produce identical results")

    assert entropy_match and heat_match, "Seeding should produce deterministic results"


def main():
    """Run all tests."""
    print("\n" + "="*70)
    print("ENTROPY EXPM_MULTIPLY VALIDATION TESTS")
    print("="*70)

    results = []

    try:
        results.append(("Small graph comparison", test_expm_multiply_vs_eigenvalue_small_graph()))
    except Exception as e:
        print(f"\n✗ Test 1 failed with exception: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Small graph comparison", False))

    try:
        results.append(("Basic functionality", test_entropy_mode_basic()))
    except Exception as e:
        print(f"\n✗ Test 2 failed with exception: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Basic functionality", False))

    try:
        results.append(("Deterministic seeding", test_deterministic_seeding()))
    except Exception as e:
        print(f"\n✗ Test 3 failed with exception: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Deterministic seeding", False))

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    for test_name, passed in results:
        status = "✓ PASS" if passed else "✗ FAIL"
        print(f"{status}: {test_name}")

    all_passed = all(passed for _, passed in results)
    if all_passed:
        print("\n✓✓✓ ALL TESTS PASSED ✓✓✓")
        return 0
    else:
        print("\n✗✗✗ SOME TESTS FAILED ✗✗✗")
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
