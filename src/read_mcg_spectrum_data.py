#!/usr/bin/env python
"""
Guide to reading MCG_SlaplSpect output data.

This script demonstrates how eigenvalues are accumulated and how to read
the different output formats from MCG_SlaplSpect.py
"""

import pickle as pk
import numpy as np
from collections import Counter
from pathlib import Path


# =============================================================================
# HOW DATA IS ACCUMULATED (Summary)
# =============================================================================

"""
MODE 1: eigval_dist
-------------------
Accumulation process:
1. Generate batch_size graphs (e.g., 50 graphs per batch)
2. Compute FULL spectrum for each graph -> arrays of eigenvalues
3. Concatenate all eigenvalues: eig_values = np.concatenate([array1, array2, ...])
   Example: 50 graphs × 30 nodes each = 1500 eigenvalues concatenated
4. Bin these eigenvalues into histogram bins (e.g., 500 bins)
5. Update Counter with bin_center -> count
6. Repeat for all batches
7. Final output: Single Counter object with {bin_center: total_count}

Storage format: Counter (dict-like)
File: mc_dist_eigval_{pflip}.pkl


MODE 2: eigvec_dist
-------------------
Accumulation process:
1. Generate batch_size graphs
2. Compute FIRST howmany eigenvectors (e.g., 5 smallest)
3. Extract components from each eigenvector
4. For eigenvector i, collect all its components across all nodes
5. Bin components for each eigenvector separately
6. Update list of Counters: [Counter_eigvec0, Counter_eigvec1, ...]
7. Final output: List of howmany Counters

Storage format: list[Counter]
File: mc_dist{howmany}_{pflip}_{backend}.pkl


MODE 3: eigvals
---------------
Accumulation process:
1. Generate graphs one by one
2. Compute FIRST howmany eigenvalues (e.g., 50 smallest)
3. Store numpy array for each graph: eigvlist.append(eigv)
4. Save checkpoint every period realizations
5. Final output: List of numpy arrays (one per graph realization)

Storage format: list[np.ndarray]
File: mc_eigvals_{pflip}_{idx}.pkl
"""


# =============================================================================
# EXAMPLE 1: Reading eigval_dist output
# =============================================================================

def read_eigval_dist(filename):
    """
    Read eigenvalue distribution (histogram).

    Parameters
    ----------
    filename : str
        Path to pickle file (e.g., 'mc_dist_eigval_0.1.pkl')

    Returns
    -------
    bin_centers : np.ndarray
        Bin center values (eigenvalues)
    counts : np.ndarray
        Number of eigenvalues in each bin
    """
    with open(filename, 'rb') as f:
        counter = pk.load(f)  # Counter object

    # Counter is a dict-like: {bin_center: count}
    bin_centers = np.array(sorted(counter.keys()))
    counts = np.array([counter[bc] for bc in bin_centers])

    return bin_centers, counts


def demo_eigval_dist():
    """Example: How to use eigval_dist output."""
    print("=" * 80)
    print("EXAMPLE 1: eigval_dist mode")
    print("=" * 80)

    # Suppose you ran:
    # python MCG_SlaplSpect.py 0.8 0.6 0.6 0.8 0.4 5 \
    #     --mode eigval_dist -p 0.1 --number_of_averages 1000

    filename = "mc_dist_eigval_0.1.pkl"

    print(f"\n1. What's in the file '{filename}':")
    print("   - A single Counter object")
    print("   - Keys: bin centers (eigenvalue values)")
    print("   - Values: counts (how many eigenvalues fell in that bin)")

    # Read the data
    try:
        bin_centers, counts = read_eigval_dist(filename)

        print(f"\n2. Data structure:")
        print(f"   - Number of bins: {len(bin_centers)}")
        print(f"   - Eigenvalue range: [{bin_centers.min():.3f}, {bin_centers.max():.3f}]")
        print(f"   - Total eigenvalues: {counts.sum()}")
        print(f"   - Average per bin: {counts.mean():.1f}")

        print(f"\n3. First 10 bins:")
        print("   Bin center | Count")
        print("   " + "-" * 30)
        for bc, c in zip(bin_centers[:10], counts[:10]):
            print(f"   {bc:10.6f} | {c:6d}")

        print("\n4. How to use this data:")
        print("   - Plot histogram: plt.plot(bin_centers, counts)")
        print("   - Normalize: density = counts / counts.sum()")
        print("   - Compare with theory")

        # Example calculations
        print("\n5. Example calculations:")
        density = counts / counts.sum()
        print(f"   - Peak density at λ = {bin_centers[np.argmax(density)]:.4f}")

        # Spectral gap (if there's a clear gap)
        threshold = 0.01 * density.max()
        gap_indices = np.where(density < threshold)[0]
        if len(gap_indices) > 0:
            print(f"   - Potential spectral gap near λ = {bin_centers[gap_indices[0]]:.4f}")

    except FileNotFoundError:
        print(f"\n   [File not found - this is just a demo]")

    print("\n" + "=" * 80)


# =============================================================================
# EXAMPLE 2: Reading eigvec_dist output
# =============================================================================

def read_eigvec_dist(filename):
    """
    Read eigenvector component distribution.

    Parameters
    ----------
    filename : str
        Path to pickle file (e.g., 'mc_dist5_0.1_scipy.pkl')

    Returns
    -------
    distributions : list[tuple]
        List of (bin_centers, counts) for each eigenvector
    """
    with open(filename, 'rb') as f:
        counters_list = pk.load(f)  # list[Counter]

    distributions = []
    for counter in counters_list:
        bin_centers = np.array(sorted(counter.keys()))
        counts = np.array([counter[bc] for bc in bin_centers])
        distributions.append((bin_centers, counts))

    return distributions


def demo_eigvec_dist():
    """Example: How to use eigvec_dist output."""
    print("=" * 80)
    print("EXAMPLE 2: eigvec_dist mode")
    print("=" * 80)

    # Suppose you ran:
    # python MCG_SlaplSpect.py 0.8 0.6 0.6 0.8 0.4 5 \
    #     --mode eigvec_dist --howmany 5 -p 0.1 --number_of_averages 500

    filename = "mc_dist5_0.1_scipy.pkl"

    print(f"\n1. What's in the file '{filename}':")
    print("   - A list of 5 Counter objects (one per eigenvector)")
    print("   - Each Counter: {component_value: count}")
    print("   - Components are from eigenvector elements across all realizations")

    try:
        distributions = read_eigvec_dist(filename)

        print(f"\n2. Data structure:")
        print(f"   - Number of eigenvectors: {len(distributions)}")

        for i, (bin_centers, counts) in enumerate(distributions):
            print(f"\n   Eigenvector {i} (e.g., {i+1}-th smallest eigenvalue):")
            print(f"      - Number of bins: {len(bin_centers)}")
            print(f"      - Component range: [{bin_centers.min():.6f}, {bin_centers.max():.6f}]")
            print(f"      - Total components: {counts.sum()}")

            # Inverse Participation Ratio (IPR) - measure of localization
            # Reconstruct approximate components from binned data
            components = []
            for bc, count in zip(bin_centers, counts):
                components.extend([bc] * int(count))
            components = np.array(components)

            ipr = np.sum(components**4) / (np.sum(components**2))**2
            print(f"      - Approx IPR: {ipr:.6f} (→1: localized, →0: delocalized)")

        print("\n3. How to use this data:")
        print("   - Plot each distribution: for bc, c in distributions[i]: plt.plot(bc, c)")
        print("   - Compare localization across eigenvectors")
        print("   - Study Anderson transition")

    except FileNotFoundError:
        print(f"\n   [File not found - this is just a demo]")

    print("\n" + "=" * 80)


# =============================================================================
# EXAMPLE 3: Reading eigvals output
# =============================================================================

def read_eigvals(filename):
    """
    Read raw eigenvalue arrays.

    Parameters
    ----------
    filename : str
        Path to pickle file (e.g., 'mc_eigvals_0.1_100.pkl')

    Returns
    -------
    eigvals_list : list[np.ndarray]
        List of eigenvalue arrays, one per graph realization
    """
    with open(filename, 'rb') as f:
        eigvals_list = pk.load(f)  # list[np.ndarray]

    return eigvals_list


def demo_eigvals():
    """Example: How to use eigvals output."""
    print("=" * 80)
    print("EXAMPLE 3: eigvals mode")
    print("=" * 80)

    # Suppose you ran:
    # python MCG_SlaplSpect.py 0.8 0.6 0.6 0.8 0.4 5 \
    #     --mode eigvals --howmany 50 -p 0.1 --number_of_averages 100

    filename = "mc_eigvals_0.1_100.pkl"

    print(f"\n1. What's in the file '{filename}':")
    print("   - A list of 100 numpy arrays")
    print("   - Each array: first 50 eigenvalues from one graph realization")
    print("   - No binning, raw data")

    try:
        eigvals_list = read_eigvals(filename)

        print(f"\n2. Data structure:")
        print(f"   - Number of realizations: {len(eigvals_list)}")
        print(f"   - Eigenvalues per realization: {len(eigvals_list[0])}")
        print(f"   - Total eigenvalues: {sum(len(ev) for ev in eigvals_list)}")

        # Stack into matrix for easier analysis
        eigvals_matrix = np.array(eigvals_list)  # shape: (100, 50)

        print(f"\n3. Example: First 3 realizations, first 5 eigenvalues:")
        print("   Realization | λ₀      λ₁      λ₂      λ₃      λ₄")
        print("   " + "-" * 60)
        for i in range(min(3, len(eigvals_list))):
            ev_str = "   ".join([f"{e:7.4f}" for e in eigvals_list[i][:5]])
            print(f"   {i:11d} | {ev_str}")

        print("\n4. Statistical analysis:")
        mean_spectrum = eigvals_matrix.mean(axis=0)
        std_spectrum = eigvals_matrix.std(axis=0)

        print(f"   Average spectrum (first 5):")
        for i in range(min(5, len(mean_spectrum))):
            print(f"      λ_{i}: {mean_spectrum[i]:.6f} ± {std_spectrum[i]:.6f}")

        print("\n5. How to use this data:")
        print("   - Compute statistics: mean, median, percentiles")
        print("   - Plot individual spectra")
        print("   - Level spacing analysis: Δ_i = λ_{i+1} - λ_i")
        print("   - Correlation analysis between eigenvalues")

        print("\n6. Example calculations:")

        # Average spectral gap
        gaps = eigvals_matrix[:, 1] - eigvals_matrix[:, 0]
        print(f"   - Average gap (λ₁ - λ₀): {gaps.mean():.6f} ± {gaps.std():.6f}")

        # Nearest-neighbor spacing distribution
        all_spacings = []
        for eigvals in eigvals_list:
            spacings = np.diff(eigvals)
            all_spacings.extend(spacings)
        all_spacings = np.array(all_spacings)

        print(f"   - Mean level spacing: {all_spacings.mean():.6f}")
        print(f"   - Spacing distribution width: {all_spacings.std():.6f}")

        # Compare to GOE/Poisson predictions
        normalized_spacings = all_spacings / all_spacings.mean()
        print(f"   - <s²>/<s>² = {(normalized_spacings**2).mean():.4f}")
        print(f"     (GOE ≈ 1.27, Poisson = 2.0)")

    except FileNotFoundError:
        print(f"\n   [File not found - this is just a demo]")

    print("\n" + "=" * 80)


# =============================================================================
# COMPLETE EXAMPLE: Full workflow
# =============================================================================

def complete_workflow_example():
    """Show complete workflow from running to analyzing."""
    print("\n" + "=" * 80)
    print("COMPLETE WORKFLOW EXAMPLE")
    print("=" * 80)

    print("""
STEP 1: Run MCG_SlaplSpect to generate data
--------------------------------------------

# Generate eigenvalue distribution
python MCG_SlaplSpect.py 0.8 0.6 0.6 0.8 0.4 6 \\
    --mode eigval_dist \\
    -p 0.1 \\
    --number_of_averages 2000 \\
    --bins_count 1000 \\
    --verbose

Output: mc_dist_eigval_0.1.pkl


STEP 2: Read and analyze the data
----------------------------------
""")

    code = '''
import pickle as pk
import numpy as np
import matplotlib.pyplot as plt

# Read the Counter object
with open('mc_dist_eigval_0.1.pkl', 'rb') as f:
    counter = pk.load(f)

# Extract histogram
bin_centers = np.array(sorted(counter.keys()))
counts = np.array([counter[bc] for bc in bin_centers])

# Normalize to get probability density
density = counts / (counts.sum() * np.diff(bin_centers).mean())

# Plot
plt.figure(figsize=(10, 6))
plt.plot(bin_centers, density, 'b-', linewidth=2, label='MC Graph')
plt.xlabel('Eigenvalue λ', fontsize=14)
plt.ylabel('Density ρ(λ)', fontsize=14)
plt.title('Spectral Density of MultiplicativeCascadeGraph (p_flip=0.1)', fontsize=16)
plt.grid(True, alpha=0.3)
plt.legend()
plt.savefig('spectral_density.pdf', dpi=300, bbox_inches='tight')
plt.show()

# Compute statistics
print(f"Total eigenvalues: {counts.sum()}")
print(f"Eigenvalue range: [{bin_centers.min():.4f}, {bin_centers.max():.4f}]")
print(f"Mean eigenvalue: {(bin_centers * counts).sum() / counts.sum():.4f}")
print(f"Peak density at: λ = {bin_centers[np.argmax(density)]:.4f}")
'''

    print(code)

    print("""

STEP 3: Advanced analysis - compare multiple pflip values
----------------------------------------------------------
""")

    advanced_code = '''
import pickle as pk
import numpy as np
import matplotlib.pyplot as plt

pflip_values = [0.0, 0.05, 0.1, 0.15, 0.2]
plt.figure(figsize=(12, 7))

for p in pflip_values:
    filename = f'mc_dist_eigval_{p:.2g}.pkl'

    with open(filename, 'rb') as f:
        counter = pk.load(f)

    bin_centers = np.array(sorted(counter.keys()))
    counts = np.array([counter[bc] for bc in bin_centers])
    density = counts / (counts.sum() * np.diff(bin_centers).mean())

    plt.plot(bin_centers, density, linewidth=2, label=f'p_flip = {p}')

plt.xlabel('Eigenvalue λ', fontsize=14)
plt.ylabel('Density ρ(λ)', fontsize=14)
plt.title('Spectral Density vs Edge Sign Disorder', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True, alpha=0.3)
plt.savefig('spectral_density_comparison.pdf', dpi=300, bbox_inches='tight')
plt.show()
'''

    print(advanced_code)
    print("\n" + "=" * 80)


# =============================================================================
# Main demonstration
# =============================================================================

def main():
    """Run all demonstrations."""
    print("\n" + "=" * 80)
    print("MCG_SlaplSpect DATA READING GUIDE")
    print("=" * 80)
    print("\nThis guide shows how eigenvalues are accumulated and how to read")
    print("the output files from MCG_SlaplSpect.py in each of the three modes.")
    print("=" * 80)

    demo_eigval_dist()
    print("\n")
    demo_eigvec_dist()
    print("\n")
    demo_eigvals()
    print("\n")
    complete_workflow_example()

    print("\n" + "=" * 80)
    print("For more information, see:")
    print("  - MCG_SLAPLSPECT_GUIDE.md")
    print("  - kernels/MCG_SlaplSpect.py (source code)")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
