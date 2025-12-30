#!/usr/bin/env python3
"""Test script for refactored MultispectralGraph hierarchy."""

import sys
sys.path.insert(0, 'lrgsglib/src')

from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph
from lrgsglib.nx_patches.Vicsek import VicsekGraph
from lrgsglib.nx_patches.DiracLattice import DiracCombGraph, DiracBrushGraph

def test_multiplicative_cascade():
    """Test MultiplicativeCascadeGraph."""
    print("Testing MultiplicativeCascadeGraph...")
    msg = MultiplicativeCascadeGraph(
        p1=0.8, p2=0.6, p3=0.6, p4=0.8,
        fraction=0.4,
        iterations=3
    )
    print(f"  {msg}")
    print(f"  Nodes: {msg.N}, Expected: {msg.get_expected_num_nodes()}")
    assert msg.N > 0, "Graph should have nodes"
    print("  ✓ MultiplicativeCascadeGraph works!\n")

def test_vicsek():
    """Test VicsekGraph."""
    print("Testing VicsekGraph...")
    msg = VicsekGraph(N=100, k=2)
    print(f"  {msg}")
    print(f"  Nodes: {msg.N}, Expected: {msg.get_expected_num_nodes()}")
    assert msg.N > 0, "Graph should have nodes"
    print("  ✓ VicsekGraph works!\n")

def test_dirac_comb():
    """Test DiracCombGraph."""
    print("Testing DiracCombGraph...")
    msg = DiracCombGraph(base_nodes=10, fiber_nodes=5)
    print(f"  {msg}")
    print(f"  Nodes: {msg.N}, Expected: {msg.get_expected_num_nodes()}")
    # Total = base_nodes + (base_nodes * fiber_nodes) = 10 + 50 = 60
    assert msg.N == 60, f"Expected 60 nodes, got {msg.N}"
    assert msg.is_dirac_structure(), "Should be a Dirac structure"
    print("  ✓ DiracCombGraph works!\n")

def test_dirac_brush():
    """Test DiracBrushGraph."""
    print("Testing DiracBrushGraph...")
    msg = DiracBrushGraph(base_x=4, base_y=4, fiber_nodes=3)
    print(f"  {msg}")
    print(f"  Nodes: {msg.N}, Expected: {msg.get_expected_num_nodes()}")
    # Total = base_nodes + (base_nodes * fiber_nodes) = 16 + 48 = 64
    assert msg.N == 64, f"Expected 64 nodes, got {msg.N}"
    assert msg.is_dirac_structure(), "Should be a Dirac structure"
    print("  ✓ DiracBrushGraph works!\n")

def test_dirac_spectrum():
    """Test Dirac spectral methods."""
    print("Testing Dirac spectral methods...")
    msg = DiracCombGraph(base_nodes=5, fiber_nodes=3)

    # Test base spectrum
    base_eigs = msg.compute_base_spectrum()
    print(f"  Base eigenvalues: {len(base_eigs)} values")
    assert len(base_eigs) == 5, f"Expected 5 base eigenvalues, got {len(base_eigs)}"

    # Test fiber Laplacian
    fiber_lap = msg.compute_fiber_laplacian()
    print(f"  Fiber Laplacian: {fiber_lap.shape}")
    assert fiber_lap.shape == (3, 3), f"Expected (3, 3), got {fiber_lap.shape}"

    # Test full spectrum
    # The separated method computes: base_nodes × fiber_nodes eigenvalues
    # This is the efficient computation exploiting the product structure
    spectrum = msg.compute_dirac_spectrum_separated()
    print(f"  Full spectrum: {len(spectrum)} values")
    expected_spectrum = 5 * 3  # base_nodes × fiber_nodes = 15
    assert len(spectrum) == expected_spectrum, f"Expected {expected_spectrum} eigenvalues, got {len(spectrum)}"

    # Test dimensions
    base_n, fiber_n = msg.get_base_fiber_dimensions()
    print(f"  Dimensions: base={base_n}, fiber={fiber_n}")
    assert base_n == 5 and fiber_n == 3, "Wrong dimensions"

    print("  ✓ Dirac spectral methods work!\n")

def main():
    """Run all tests."""
    print("=" * 60)
    print("Testing Refactored MultispectralGraph Hierarchy")
    print("=" * 60 + "\n")

    try:
        test_multiplicative_cascade()
        test_vicsek()
        test_dirac_comb()
        test_dirac_brush()
        test_dirac_spectrum()

        print("=" * 60)
        print("✓ ALL TESTS PASSED!")
        print("=" * 60)
        return 0

    except Exception as e:
        print("\n" + "=" * 60)
        print(f"✗ TEST FAILED: {e}")
        print("=" * 60)
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
