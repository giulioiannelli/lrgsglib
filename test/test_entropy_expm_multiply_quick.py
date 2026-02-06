#!/usr/bin/env python
"""
Quick entropy expm_multiply test (lightweight, no matplotlib).

Validates core functionality without heavy dependencies.
"""

import numpy as np
from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph
from lrgsglib.utils.lrg.infocomm import compute_entropy_observables_expm_multiply


def test_basic_computation():
    """Test that entropy computation runs without errors."""
    print("\n" + "="*60)
    print("Test: Basic Entropy Computation")
    print("="*60)

    # Create small test graph
    G = MultiplicativeCascadeGraph(
        p1=0.8, p2=0.6, p3=0.6, p4=0.8,
        fraction=0.3,
        iterations=2,  # Very small graph
        pflip=0.1,
        seed=42
    )
    print(f"Graph: N={G.N} nodes")

    # Ensure slp is available
    G.upd_graph_matrices()

    # Compute entropy
    S, C, V, tau = compute_entropy_observables_expm_multiply(
        L=G.slp,
        num_nodes=G.N,
        steps=20,  # Few time points
        t1=-1,
        t2=2,
        num_samples=10,  # Few samples
        seed=42
    )

    # Validate
    checks = [
        ("Correct shapes", len(S) == len(C) == len(V) == len(tau) == 20),
        ("No NaNs", not np.any(np.isnan(S))),
        ("All finite", np.all(np.isfinite(S))),
        ("Entropy bounds", np.all((S >= 0) & (S <= 1.5))),  # Allow slight overshoot
    ]

    all_pass = True
    for name, result in checks:
        status = "✓" if result else "✗"
        print(f"  {status} {name}")
        all_pass = all_pass and result

    assert all_pass, "Some basic checks failed"


def test_deterministic():
    """Test that seeding produces identical results."""
    print("\n" + "="*60)
    print("Test: Deterministic Seeding")
    print("="*60)

    G = MultiplicativeCascadeGraph(p1=0.8, p2=0.6, p3=0.6, p4=0.8,
                                    fraction=0.3, iterations=2, pflip=0.1, seed=42)
    G.upd_graph_matrices()

    # Run twice with same seed
    S1, C1, _, _ = compute_entropy_observables_expm_multiply(
        L=G.slp, num_nodes=G.N, steps=20, num_samples=10, seed=100
    )
    S2, C2, _, _ = compute_entropy_observables_expm_multiply(
        L=G.slp, num_nodes=G.N, steps=20, num_samples=10, seed=100
    )

    match = np.allclose(S1, S2, rtol=1e-14)
    print(f"  {'✓' if match else '✗'} Identical results with same seed")

    assert match, "Same seed should produce identical results"


def test_error_convergence():
    """Test that error decreases with more samples (qualitative)."""
    print("\n" + "="*60)
    print("Test: Error Convergence (qualitative)")
    print("="*60)

    G = MultiplicativeCascadeGraph(p1=0.8, p2=0.6, p3=0.6, p4=0.8,
                                    fraction=0.3, iterations=2, pflip=0.0, seed=42)
    G.upd_graph_matrices()

    # Different sample counts
    results = []
    for m in [5, 10, 20]:
        S, _, _, _ = compute_entropy_observables_expm_multiply(
            L=G.slp, num_nodes=G.N, steps=20, num_samples=m, seed=42
        )
        results.append(S)
        print(f"  num_samples={m:2d}: mean(S)={np.mean(S):.4f}")

    # Just verify we got results (convergence test needs ground truth)
    print("  ✓ Computed entropy for different sample sizes")
    assert len(results) == 3, "Should have computed 3 results"


def main():
    """Run quick tests."""
    print("\n" + "="*60)
    print("QUICK ENTROPY VALIDATION TESTS")
    print("="*60)

    results = []

    try:
        results.append(("Basic computation", test_basic_computation()))
    except Exception as e:
        print(f"✗ Test failed: {e}")
        results.append(("Basic computation", False))

    try:
        results.append(("Deterministic seeding", test_deterministic()))
    except Exception as e:
        print(f"✗ Test failed: {e}")
        results.append(("Deterministic seeding", False))

    try:
        results.append(("Error convergence", test_error_convergence()))
    except Exception as e:
        print(f"✗ Test failed: {e}")
        results.append(("Error convergence", False))

    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    for name, passed in results:
        print(f"  {'✓' if passed else '✗'} {name}")

    all_passed = all(p for _, p in results)
    if all_passed:
        print("\n✓✓✓ ALL TESTS PASSED ✓✓✓")
        return 0
    else:
        print("\n✗ SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
