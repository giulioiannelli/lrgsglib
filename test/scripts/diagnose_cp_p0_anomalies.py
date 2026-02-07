#!/usr/bin/env python3
"""
Diagnostic script to find p=0.0 contact process runs that don't decay to zero.
Focus on gamma values around 0.5 that show anomalous behavior.

This script checks contact process simulation data files for runs that should
decay to zero but don't (particularly in the subcritical and critical regimes).
"""
from pathlib import Path
import numpy as np
import re
import os
import argparse
from typing import List, Tuple

# Default configuration
DEFAULT_DATA_DIR = Path(__file__).parent.parent.parent / "data" / "cluster_data" / "new_cp_filetransfer" / "l2d_squared" / "cntct" / "N=9216"
PFLIP = "0"  # Use string to match filenames exactly
THRESHOLD = 1e-6  # Density threshold to consider as "not zero"
NS = 500  # Number of samples per run

def extract_gamma_from_filename(filepath: Path) -> float:
    """Extract gamma value from filename."""
    match = re.search(r"gamma=([^_]+)", filepath.name)
    if match:
        return float(match.group(1))
    raise ValueError(f"Could not extract gamma from {filepath.name}")

def extract_na_from_filename(filepath: Path) -> int:
    """Extract number of runs (na) from filename."""
    match = re.search(r"na=(\d+)", filepath.name)
    if match:
        return int(match.group(1))
    raise ValueError(f"Could not extract na from {filepath.name}")

def check_file(filepath: Path, threshold: float = THRESHOLD) -> Tuple[float, int, float, float, int]:
    """
    Check if runs in this file decay to zero.

    Args:
        filepath: Path to the binary data file
        threshold: Density threshold to consider as "not zero"

    Returns:
        gamma: gamma value
        na: number of runs
        mean_final: mean final density across all runs
        max_final: maximum final density across all runs
        num_nonzero: number of runs with non-zero final density
    """
    gamma = extract_gamma_from_filename(filepath)
    na = extract_na_from_filename(filepath)

    try:
        # Try to load with expected shape
        densitydata = np.fromfile(filepath, dtype=np.float64).reshape(na, -1)
    except ValueError:
        # If that fails, try reverse
        densitydata = np.fromfile(filepath, dtype=np.float64).reshape(-1, NS)
        na = densitydata.shape[0]

    # Get final densities for all runs
    final_densities = densitydata[:, -1]

    mean_final = final_densities.mean()
    max_final = final_densities.max()
    num_nonzero = np.sum(final_densities > threshold)

    return gamma, na, mean_final, max_final, num_nonzero

def main():
    parser = argparse.ArgumentParser(
        description="Diagnose contact process p=0.0 runs that don't decay to zero"
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help=f"Path to data directory (default: {DEFAULT_DATA_DIR})",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=THRESHOLD,
        help=f"Density threshold to consider as non-zero (default: {THRESHOLD})",
    )
    args = parser.parse_args()

    data_dir = args.data_dir
    threshold = args.threshold

    # Check if directory exists
    if not data_dir.exists():
        print(f"❌ Error: Data directory not found: {data_dir}")
        print(f"   Please specify the correct path with --data-dir")
        return 1

    # Find all p=0 files
    pattern = f"dens_p={PFLIP}_gamma=*_uniform_rand_*.bin"
    files = sorted(data_dir.glob(pattern))

    print(f"Data directory: {data_dir}")
    print(f"Found {len(files)} files for p={PFLIP}")
    print(f"Checking for runs that don't decay below {threshold}...")
    print("=" * 80)

    # Collect results
    anomalies = []
    all_results = []

    for filepath in files:
        gamma, na, mean_final, max_final, num_nonzero = check_file(filepath, threshold)
        all_results.append((gamma, na, mean_final, max_final, num_nonzero, filepath))

        if mean_final > threshold or num_nonzero > 0:
            anomalies.append((gamma, na, mean_final, max_final, num_nonzero, filepath))

    # Sort results by gamma
    all_results.sort(key=lambda x: x[0])
    anomalies.sort(key=lambda x: x[0])

    # Display anomalies
    if anomalies:
        print(f"\n🚨 ANOMALIES FOUND: {len(anomalies)} files with non-zero final density")
        print("=" * 80)
        for gamma, na, mean_final, max_final, num_nonzero, filepath in anomalies:
            print(f"\nγ = {gamma:.4f}")
            print(f"  File: {filepath.name}")
            print(f"  Number of runs: {na}")
            print(f"  Mean final density: {mean_final:.6e}")
            print(f"  Max final density: {max_final:.6e}")
            print(f"  Runs with non-zero final: {num_nonzero}/{na} ({100*num_nonzero/na:.1f}%)")

            # Load and show which specific runs are problematic
            densitydata = np.fromfile(filepath, dtype=np.float64)
            try:
                densitydata = densitydata.reshape(na, -1)
            except ValueError:
                densitydata = densitydata.reshape(-1, NS)

            final_densities = densitydata[:, -1]
            nonzero_indices = np.where(final_densities > threshold)[0]

            if len(nonzero_indices) <= 10:
                print(f"  Problematic run indices: {list(nonzero_indices)}")
                print(f"  Their final densities: {final_densities[nonzero_indices]}")
            else:
                print(f"  Problematic run indices (first 10): {list(nonzero_indices[:10])}")
                print(f"  Their final densities (first 10): {final_densities[nonzero_indices[:10]]}")
    else:
        print("\n✓ No anomalies found - all runs decay to zero properly")

    # Show gamma range around 0.5
    print("\n" + "=" * 80)
    print("GAMMA VALUES AROUND 0.5:")
    print("=" * 80)
    for gamma, na, mean_final, max_final, num_nonzero, filepath in all_results:
        if 0.49 <= gamma <= 0.53:
            status = "🚨 ANOMALY" if mean_final > threshold else "✓ OK"
            print(f"{status}  γ={gamma:.4f}  mean_final={mean_final:.2e}  max_final={max_final:.2e}  nonzero={num_nonzero}/{na}")

    # Summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY:")
    print("=" * 80)
    gamma_values = [r[0] for r in all_results]
    if gamma_values:
        print(f"Gamma range: {min(gamma_values):.4f} to {max(gamma_values):.4f}")
    print(f"Total files analyzed: {len(all_results)}")
    print(f"Files with anomalies: {len(anomalies)}")

    if anomalies:
        anomaly_gammas = [a[0] for a in anomalies]
        print(f"Anomaly gamma range: {min(anomaly_gammas):.4f} to {max(anomaly_gammas):.4f}")

if __name__ == "__main__":
    main()
