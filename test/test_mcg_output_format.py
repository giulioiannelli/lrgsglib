"""
Test script to verify MCG_SlaplSpect output file format.

Verifies that saved eigenvalue files:
1. Have correct structure (list[np.ndarray])
2. Can be loaded and processed
3. Work with entropy calculation functions
"""

import subprocess
import pickle
import sys
from pathlib import Path
import numpy as np
import pytest

ROOT_DIR = Path(__file__).resolve().parents[1]
MCG_SCRIPT = ROOT_DIR / "src" / "MCG_SlaplSpect.py"

# Add src to path
sys.path.insert(0, str(ROOT_DIR / "src"))

from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph
from lrgsglib.utils.lrg.infocomm import compute_entropy_observables_from_eigenvalues

pytestmark = [pytest.mark.integration, pytest.mark.slow]

def get_output_dir(fraction: float, iterations: int) -> Path:
    """Get the output directory for test files."""
    G = MultiplicativeCascadeGraph(
        p1=0.8, p2=0.6, p3=0.6, p4=0.8,
        fraction=fraction,
        iterations=iterations,
        only_const_mode=True
    )
    return G.path_spect


def test_file_format():
    """Test that output files have correct format and work with entropy functions."""
    print("\n" + "="*70)
    print("Testing MCG_SlaplSpect Output File Format")
    print("="*70)

    # First, run the program to generate a test file
    print("\nStep 1: Generating test eigenvalue file...")
    cmd = [
        "python", str(MCG_SCRIPT),
        "0.8", "0.6", "0.6", "0.8",  # p1, p2, p3, p4
        "0.7",                        # fraction
        "5",                          # iterations (716 nodes with frac=0.7)
        "--pflip", "0.05",            # pflip
        "-na", "4",                   # number_of_averages
        "-m", "eigvals",              # mode
        "--backend", "scipy",         # backend
        "--howmany", "0",             # full spectrum
        "--save_frequency", "2",       # save every 2 iterations
        "--verbose"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path.cwd())

    if result.returncode != 0:
        print(f"❌ FAILED: Program execution failed")
        print(f"STDERR: {result.stderr}")
        pytest.fail("Program execution failed")

    print("✓ Program ran successfully")

    # Load the saved file
    output_dir = get_output_dir(fraction=0.7, iterations=5)
    test_file = output_dir / "mc_eigvals_p=0.05_na=4.pkl"  # pflip=0.05, navg=4

    if not test_file.exists():
        print(f"❌ FAILED: Output file not found: {test_file}")
        print("Available files:")
        for f in output_dir.glob("mc_eigvals_*"):
            print(f"  - {f.name}")
        pytest.fail(f"Output file not found: {test_file}")

    print(f"\nStep 2: Loading eigenvalue file: {test_file.name}")

    try:
        with open(test_file, 'rb') as f:
            eigvlist = pickle.load(f)
    except Exception as e:
        print(f"❌ FAILED: Could not load file: {e}")
        pytest.fail(f"Could not load file: {e}")

    print("✓ File loaded successfully")

    # Verify structure
    print("\nStep 3: Verifying file structure...")

    if not isinstance(eigvlist, list):
        print(f"❌ FAILED: Expected list, got {type(eigvlist)}")
        assert False, f"Expected list, got {type(eigvlist)}"

    print(f"✓ File contains a list with {len(eigvlist)} elements")

    if len(eigvlist) != 4:
        print(f"❌ FAILED: Expected 4 realizations, got {len(eigvlist)}")
        assert False, f"Expected 4 realizations, got {len(eigvlist)}"

    print(f"✓ Correct number of realizations: {len(eigvlist)}")

    for i, eigv in enumerate(eigvlist):
        if not isinstance(eigv, np.ndarray):
            print(f"❌ FAILED: Element {i} is not np.ndarray: {type(eigv)}")
            assert False, f"Element {i} is not np.ndarray: {type(eigv)}"

        if eigv.ndim != 1:
            print(f"❌ FAILED: Element {i} has wrong dimensions: {eigv.ndim} (expected 1)")
            assert False, f"Element {i} has wrong dimensions: {eigv.ndim}"

    print(f"✓ All elements are 1D numpy arrays")

    # Print details
    print("\nEigenvalue array details:")
    for i, eigv in enumerate(eigvlist):
        print(f"  Realization {i}: shape={eigv.shape}, min={eigv.min():.6f}, max={eigv.max():.6f}, mean={eigv.mean():.6f}")

    # Test entropy computation
    print("\nStep 4: Testing entropy computation with saved eigenvalues...")

    try:
        for i, eigvals in enumerate(eigvlist):
            N = len(eigvals)

            entropy, deriv, variance, tau = compute_entropy_observables_from_eigenvalues(
                eigenvalues=eigvals,
                num_nodes=N,
                steps=50,
                t1=-2,
                t2=3,
                typf=np.float64,
                threshold=1e-10
            )

            # Verify outputs
            if len(entropy) != 50:
                print(f"❌ FAILED: Entropy array has wrong length: {len(entropy)} (expected 50)")
                assert False, f"Entropy array has wrong length: {len(entropy)}"

            if len(deriv) != 50:
                print(f"❌ FAILED: Derivative array has wrong length: {len(deriv)} (expected 50)")
                assert False, f"Derivative array has wrong length: {len(deriv)}"

            if len(tau) != 50:
                print(f"❌ FAILED: Time grid has wrong length: {len(tau)} (expected 50)")
                assert False, f"Time grid has wrong length: {len(tau)}"

            # Check for valid values (no NaN/Inf)
            if np.any(np.isnan(entropy)) or np.any(np.isinf(entropy)):
                print(f"❌ FAILED: Entropy contains NaN or Inf values")
                assert False, "Entropy contains NaN or Inf values"

            print(f"  Realization {i}: entropy range=[{entropy.min():.4f}, {entropy.max():.4f}], " +
                  f"specific_heat range=[{deriv.min():.4f}, {deriv.max():.4f}]")

    except Exception as e:
        print(f"❌ FAILED: Entropy computation failed: {e}")
        import traceback
        traceback.print_exc()
        pytest.fail(f"Entropy computation failed: {e}")

    print("✓ All entropy computations succeeded")

    # Test averaging across realizations
    print("\nStep 5: Testing averaging across realizations...")

    try:
        all_entropies = []
        all_specific_heats = []

        for eigvals in eigvlist:
            N = len(eigvals)
            entropy, deriv, variance, tau = compute_entropy_observables_from_eigenvalues(
                eigenvalues=eigvals,
                num_nodes=N,
                steps=50,
                t1=-2,
                t2=3,
                typf=np.float64,
                threshold=1e-10
            )
            all_entropies.append(entropy)
            all_specific_heats.append(deriv)

        # Average
        mean_entropy = np.mean(all_entropies, axis=0)
        std_entropy = np.std(all_entropies, axis=0)
        mean_specific_heat = np.mean(all_specific_heats, axis=0)
        std_specific_heat = np.std(all_specific_heats, axis=0)

        print(f"  Mean entropy range: [{mean_entropy.min():.4f}, {mean_entropy.max():.4f}]")
        print(f"  Mean entropy std: [{std_entropy.min():.4f}, {std_entropy.max():.4f}]")
        print(f"  Mean specific heat range: [{mean_specific_heat.min():.4f}, {mean_specific_heat.max():.4f}]")
        print(f"  Mean specific heat std: [{std_specific_heat.min():.4f}, {std_specific_heat.max():.4f}]")

    except Exception as e:
        print(f"❌ FAILED: Averaging failed: {e}")
        pytest.fail(f"Averaging failed: {e}")

    print("✓ Averaging across realizations works correctly")

    # Clean up
    print("\nCleaning up test file...")
    test_file.unlink()
    print(f"✓ Removed: {test_file}")

    print("\n" + "="*70)
    print("✅ FILE FORMAT TEST PASSED")
    print("="*70)
    print("\nSummary:")
    print("  ✓ File has correct structure: list[np.ndarray]")
    print("  ✓ Eigenvalues can be loaded without errors")
    print("  ✓ Entropy observables can be computed from eigenvalues")
    print("  ✓ Averaging across realizations works")
    print("\nConclusion: Output files are in the correct format for entropy analysis!")

    return None


def main():
    """Run the file format test."""
    print("\n" + "#"*70)
    print("# MCG_SlaplSpect Output File Format Verification")
    print("#"*70)

    try:
        test_file_format()
    except AssertionError:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
