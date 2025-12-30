#!/usr/bin/env python3
"""
Test script to compare topology-based vs edge-sign-based clustering
for SCS order parameter calculations.

This runs sweeps for N=[8,16,32,64] with J0 from 0.9 to 2.0 to check
if the std peak location shifts with system size.
"""

import sys
from pathlib import Path

# Add lrgsglib/src to path
repo_root = Path("..") / "lrgsglib"
sys.path.insert(0, str(repo_root / "src"))

import numpy as np
import matplotlib.pyplot as plt
from SCS_TransCluster_sweep import run_batch

# Test parameters
N_list = [8, 16, 32, 64]
g_list = [2.0]  # g=2 for comparison
j0_values = np.linspace(0.9, 2.0, 50)  # Focus on transition region

gamma = 1.0
J = 1.0
diagonal = -1.0
navg = 2000  # Production run with higher statistics

backend = "numpy"
float_type = "float64"
partition_rule = ">0"
nproc = 10
seed = None  # No fixed seed - different realizations for statistical averaging
verbose = True

print("=" * 70)
print("Testing Edge-Sign-Based vs Topology-Based Clustering")
print("=" * 70)
print(f"\nParameters:")
print(f"  N values: {N_list}")
print(f"  g values: {g_list}")
print(f"  J0 range: {j0_values[0]:.2f} to {j0_values[-1]:.2f} ({len(j0_values)} points)")
print(f"  navg: {navg}")
print(f"  gamma: {gamma}, J: {J}, diagonal: {diagonal}")
print()

# Run topology-based (original method)
print("\n" + "=" * 70)
print("1. Running TOPOLOGY-BASED clustering (original method)")
print("=" * 70)
workdir_topo = "test_scs_topology_clustering"
run_batch(
    N_list=N_list,
    g_list=g_list,
    j0_values=j0_values,
    gamma=gamma,
    J=J,
    diagonal=diagonal,
    navg=navg,
    backend=backend,
    float_type=float_type,
    partition_rule=partition_rule,
    workdir=workdir_topo,
    out_suffix=None,
    nproc=nproc,
    seed=seed,
    verbose=verbose,
    progress=True,
    use_edge_sign_clustering=False,  # Original method
)

# Run edge-sign-based (new method)
print("\n" + "=" * 70)
print("2. Running EDGE-SIGN-BASED clustering (new method)")
print("=" * 70)
workdir_edge = "test_scs_edgesign_clustering"
run_batch(
    N_list=N_list,
    g_list=g_list,
    j0_values=j0_values,
    gamma=gamma,
    J=J,
    diagonal=diagonal,
    navg=navg,
    backend=backend,
    float_type=float_type,
    partition_rule=partition_rule,
    workdir=workdir_edge,
    out_suffix=None,
    nproc=nproc,
    seed=seed,
    verbose=verbose,
    progress=True,
    use_edge_sign_clustering=True,  # New method
)

print("\n" + "=" * 70)
print("Sweep complete! Run the analysis script to visualize results.")
print("=" * 70)
