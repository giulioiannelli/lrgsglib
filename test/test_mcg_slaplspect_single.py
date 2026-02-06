"""
Test script for MCG_SlaplSpect single program execution.

Tests all three computation modes with small parameters to verify:
1. Programs run without errors
2. Output files are created
3. Files can be loaded and have correct structure
"""

import subprocess
import pickle
import os
from pathlib import Path
import sys
import pytest

ROOT_DIR = Path(__file__).resolve().parents[1]
MCG_SCRIPT = ROOT_DIR / "src" / "MCG_SlaplSpect.py"

# Add src to path for imports
sys.path.insert(0, str(ROOT_DIR / "src"))

from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph

pytestmark = [pytest.mark.integration, pytest.mark.slow]


def get_output_dir():
    """Get the output directory for test files."""
    # Create a temporary test MCG to get path_spect
    G = MultiplicativeCascadeGraph(
        p1=0.8, p2=0.6, p3=0.6, p4=0.8,
        fraction=0.5,
        iterations=5,
        only_const_mode=True
    )
    return G.path_spect


def test_eigvals_mode():
    """Test raw eigenvalue mode with small system."""
    print("\n" + "="*70)
    print("Testing eigvals mode (raw eigenvalues for entropy calculations)")
    print("="*70)

    cmd = [
        "python", str(MCG_SCRIPT),
        "0.8", "0.6", "0.6", "0.8",  # p1, p2, p3, p4
        "0.5",                        # fraction
        "5",                          # iterations (small: 512 nodes)
        "--pflip", "0.1",             # pflip
        "-na", "3",                   # number_of_averages (small for test)
        "-m", "eigvals",              # mode
        "--backend", "scipy",         # backend
        "--howmany", "0",             # full spectrum
        "--save_frequency", "1",      # save every iteration
        "--verbose"
    ]

    print(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path.cwd())

    print(f"\nReturn code: {result.returncode}")
    if result.stdout:
        print(f"STDOUT:\n{result.stdout}")
    if result.stderr:
        print(f"STDERR:\n{result.stderr}")

    if result.returncode != 0:
        print(f"❌ FAILED: eigvals mode returned non-zero exit code")
        pytest.fail("eigvals mode returned non-zero exit code")

    # Verify output file exists
    output_dir = get_output_dir()
    expected_file = output_dir / "mc_eigvals_p=0.1_na=3.pkl"  # pflip=0.1, navg=3

    if not expected_file.exists():
        print(f"❌ FAILED: Expected output file not found: {expected_file}")
        pytest.fail(f"Expected output file not found: {expected_file}")

    # Load and verify structure
    try:
        with open(expected_file, 'rb') as f:
            eigvlist = pickle.load(f)

        print(f"\n✓ Output file loaded successfully")
        print(f"  - Type: {type(eigvlist)}")
        print(f"  - Length: {len(eigvlist)} realizations")
        print(f"  - First eigenvalue array shape: {eigvlist[0].shape}")
        print(f"  - First eigenvalue array type: {type(eigvlist[0])}")

        # Verify structure
        if not isinstance(eigvlist, list):
            print(f"❌ FAILED: Expected list, got {type(eigvlist)}")
            assert False, f"Expected list, got {type(eigvlist)}"

        if len(eigvlist) != 3:
            print(f"❌ FAILED: Expected 3 realizations, got {len(eigvlist)}")
            assert False, f"Expected 3 realizations, got {len(eigvlist)}"

        import numpy as np
        for i, eigv in enumerate(eigvlist):
            if not isinstance(eigv, np.ndarray):
                print(f"❌ FAILED: Realization {i} is not np.ndarray: {type(eigv)}")
                assert False, f"Realization {i} is not np.ndarray: {type(eigv)}"

        print(f"✓ File structure is correct: list[np.ndarray]")

        # Clean up test file
        expected_file.unlink()
        print(f"✓ Cleaned up test file: {expected_file}")

        print("\n✅ eigvals mode test PASSED")
        return None

    except Exception as e:
        print(f"❌ FAILED: Error loading/verifying file: {e}")
        pytest.fail(f"Error loading/verifying file: {e}")


def test_eigval_dist_mode():
    """Test eigenvalue distribution mode."""
    print("\n" + "="*70)
    print("Testing eigval_dist mode (eigenvalue histogram)")
    print("="*70)

    cmd = [
        "python", str(MCG_SCRIPT),
        "0.8", "0.6", "0.6", "0.8",  # p1, p2, p3, p4
        "0.5",                        # fraction
        "5",                          # iterations
        "--pflip", "0.0",             # pflip
        "-na", "5",                   # number_of_averages
        "-m", "eigval_dist",          # mode
        "--backend", "scipy",         # backend
        "-bc", "50",                  # bins_count
        "--verbose"
    ]

    print(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path.cwd())

    print(f"\nReturn code: {result.returncode}")
    if result.stdout:
        print(f"STDOUT:\n{result.stdout}")
    if result.stderr:
        print(f"STDERR:\n{result.stderr}")

    if result.returncode != 0:
        print(f"❌ FAILED: eigval_dist mode returned non-zero exit code")
        pytest.fail("eigval_dist mode returned non-zero exit code")

    # Verify output file
    output_dir = get_output_dir()
    expected_file = output_dir / "mc_dist_eigval_p=0_na=5.pkl"  # pflip=0.0, navg=5

    if not expected_file.exists():
        print(f"❌ FAILED: Expected output file not found: {expected_file}")
        # List files in directory to help debug
        print(f"Files in {output_dir}:")
        for f in output_dir.glob("mc_dist_eigval_*"):
            print(f"  - {f.name}")
        pytest.fail(f"Expected output file not found: {expected_file}")

    # Load and verify
    try:
        with open(expected_file, 'rb') as f:
            data = pickle.load(f)

        print(f"\n✓ Output file loaded successfully")
        print(f"  - Type: {type(data)}")

        # Clean up
        expected_file.unlink()
        print(f"✓ Cleaned up test file: {expected_file}")

        print("\n✅ eigval_dist mode test PASSED")
        return None

    except Exception as e:
        print(f"❌ FAILED: Error loading file: {e}")
        pytest.fail(f"Error loading file: {e}")


def test_eigvec_dist_mode():
    """Test eigenvector distribution mode."""
    print("\n" + "="*70)
    print("Testing eigvec_dist mode (eigenvector component histogram)")
    print("="*70)

    cmd = [
        "python", str(MCG_SCRIPT),
        "0.8", "0.6", "0.6", "0.8",  # p1, p2, p3, p4
        "0.5",                        # fraction
        "5",                          # iterations
        "--pflip", "0.0",             # pflip
        "-na", "5",                   # number_of_averages
        "-m", "eigvec_dist",          # mode
        "--backend", "scipy",         # backend
        "-em", "smallest",            # eigen_mode (selection mode)
        "-hm", "2",                   # howmany (first 2 eigenvectors)
        "-bc", "50",                  # bins_count
        "--verbose"
    ]

    print(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=Path.cwd())

    print(f"\nReturn code: {result.returncode}")
    if result.stdout:
        print(f"STDOUT:\n{result.stdout}")
    if result.stderr:
        print(f"STDERR:\n{result.stderr}")

    if result.returncode != 0:
        print(f"❌ FAILED: eigvec_dist mode returned non-zero exit code")
        pytest.fail("eigvec_dist mode returned non-zero exit code")

    # Verify output file
    output_dir = get_output_dir()
    expected_file = output_dir / "mc_dist2_p=0_smallest_na=5.pkl"  # howmany=2, pflip=0.0, mode=smallest, navg=5

    if not expected_file.exists():
        print(f"❌ FAILED: Expected output file not found: {expected_file}")
        # List files in directory to help debug
        print(f"Files in {output_dir}:")
        for f in output_dir.glob("mc_dist*"):
            print(f"  - {f.name}")
        pytest.fail(f"Expected output file not found: {expected_file}")

    # Load and verify
    try:
        with open(expected_file, 'rb') as f:
            data = pickle.load(f)

        print(f"\n✓ Output file loaded successfully")
        print(f"  - Type: {type(data)}")

        # Clean up
        expected_file.unlink()
        print(f"✓ Cleaned up test file: {expected_file}")

        print("\n✅ eigvec_dist mode test PASSED")
        return None

    except Exception as e:
        print(f"❌ FAILED: Error loading file: {e}")
        pytest.fail(f"Error loading file: {e}")


def main():
    """Run all tests."""
    print("\n" + "#"*70)
    print("# MCG_SlaplSpect Single Program Tests")
    print("#"*70)

    try:
        test_eigvals_mode()
        test_eigval_dist_mode()
        test_eigvec_dist_mode()
    except AssertionError:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
