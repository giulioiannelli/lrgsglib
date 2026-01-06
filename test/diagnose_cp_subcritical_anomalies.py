#!/usr/bin/env python3
"""
Find anomalous runs in the SUBCRITICAL regime (gamma < 0.5) that don't decay to zero.

This diagnostic script analyzes contact process simulation data in the subcritical
regime where ALL runs should decay to zero (for p=0.0). Any runs that don't decay
represent a serious problem with the simulation or data.
"""
from pathlib import Path
import numpy as np
import re
import argparse
from typing import Dict, List, Any

# Default configuration
DEFAULT_DATA_DIR = Path(__file__).parent.parent.parent / "data" / "cluster_data" / "new_cp_filetransfer" / "l2d_squared" / "cntct" / "N=9216"
PFLIP = "0"
THRESHOLD = 1e-6
GAMMA_CRITICAL = 0.5  # Critical point for contact process on square lattice
NS = 500

def extract_gamma_from_filename(filepath: Path) -> float:
    match = re.search(r"gamma=([^_]+)", filepath.name)
    if match:
        return float(match.group(1))
    raise ValueError(f"Could not extract gamma from {filepath.name}")

def extract_na_from_filename(filepath: Path) -> int:
    match = re.search(r"na=(\d+)", filepath.name)
    if match:
        return int(match.group(1))
    raise ValueError(f"Could not extract na from {filepath.name}")

def main():
    parser = argparse.ArgumentParser(
        description="Diagnose contact process subcritical regime anomalies (γ < 0.5, p=0.0)"
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
    parser.add_argument(
        "--gamma-critical",
        type=float,
        default=GAMMA_CRITICAL,
        help=f"Critical gamma value (default: {GAMMA_CRITICAL})",
    )
    args = parser.parse_args()

    data_dir = args.data_dir
    threshold = args.threshold
    gamma_critical = args.gamma_critical

    # Check if directory exists
    if not data_dir.exists():
        print(f"❌ Error: Data directory not found: {data_dir}")
        print(f"   Please specify the correct path with --data-dir")
        return 1

    pattern = f"dens_p={PFLIP}_gamma=*_uniform_rand_*.bin"
    files = sorted(data_dir.glob(pattern))

    print("=" * 80)
    print("SUBCRITICAL REGIME ANOMALY ANALYSIS (γ < 0.5)")
    print("=" * 80)
    print(f"Data directory: {data_dir}")
    print(f"For p=0.0, ALL runs should decay to zero when γ < {gamma_critical}")
    print(f"Checking for runs with final density > {threshold}...")
    print("=" * 80)

    subcritical_anomalies: List[Dict[str, Any]] = []

    for filepath in files:
        gamma = extract_gamma_from_filename(filepath)

        # Only analyze subcritical regime
        if gamma >= gamma_critical:
            continue

        na = extract_na_from_filename(filepath)

        try:
            densitydata = np.fromfile(filepath, dtype=np.float64).reshape(na, -1)
        except ValueError:
            densitydata = np.fromfile(filepath, dtype=np.float64).reshape(-1, NS)
            na = densitydata.shape[0]

        final_densities = densitydata[:, -1]
        nonzero_indices = np.where(final_densities > threshold)[0]

        if len(nonzero_indices) > 0:
            subcritical_anomalies.append({
                'gamma': gamma,
                'filepath': filepath,
                'na': na,
                'anomalous_runs': nonzero_indices,
                'final_densities': final_densities[nonzero_indices],
                'densitydata': densitydata,
            })

    if not subcritical_anomalies:
        print("\n✓ No anomalies found in subcritical regime!")
        return

    print(f"\n🚨 FOUND {len(subcritical_anomalies)} FILES WITH ANOMALIES IN SUBCRITICAL REGIME\n")

    for i, anom in enumerate(subcritical_anomalies, 1):
        print(f"\n{'='*80}")
        print(f"ANOMALY {i}: γ = {anom['gamma']:.4f} (SUBCRITICAL)")
        print(f"{'='*80}")
        print(f"File: {anom['filepath'].name}")
        print(f"Total runs: {anom['na']}")
        print(f"Anomalous runs: {len(anom['anomalous_runs'])}/{anom['na']} ({100*len(anom['anomalous_runs'])/anom['na']:.2f}%)")
        print(f"Anomalous run indices: {list(anom['anomalous_runs'])}")
        print(f"Their final densities: {anom['final_densities']}")

        # Check if these runs stay at 1.0 throughout or if they evolve
        for run_idx in anom['anomalous_runs'][:3]:  # Check first 3
            trajectory = anom['densitydata'][run_idx, :]
            print(f"\n  Run {run_idx} trajectory:")
            print(f"    Start: {trajectory[0]:.6f}")
            print(f"    t=100: {trajectory[100]:.6f}")
            print(f"    t=200: {trajectory[200]:.6f}")
            print(f"    t=300: {trajectory[300]:.6f}")
            print(f"    t=400: {trajectory[400]:.6f}")
            print(f"    Final: {trajectory[-1]:.6f}")

            # Check if it's stuck at 1.0
            if np.all(np.isclose(trajectory, 1.0)):
                print(f"    ⚠️  STUCK AT 1.0 FOR ENTIRE RUN!")
            elif np.all(trajectory[-50:] > 0.99):
                print(f"    ⚠️  SATURATED AT 1.0 (shouldn't happen in subcritical!)")
            else:
                print(f"    Density range: [{trajectory.min():.6f}, {trajectory.max():.6f}]")

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY OF SUBCRITICAL ANOMALIES:")
    print("=" * 80)
    total_anomalous_runs = sum(len(a['anomalous_runs']) for a in subcritical_anomalies)
    total_runs = sum(a['na'] for a in subcritical_anomalies)
    print(f"Total files with anomalies: {len(subcritical_anomalies)}")
    print(f"Total anomalous runs: {total_anomalous_runs}/{total_runs} ({100*total_anomalous_runs/total_runs:.2f}%)")
    print(f"\nFiles:")
    for anom in subcritical_anomalies:
        print(f"  - {anom['filepath'].name}")
        print(f"    γ={anom['gamma']:.4f}, {len(anom['anomalous_runs'])}/{anom['na']} anomalous runs")

if __name__ == "__main__":
    main()
