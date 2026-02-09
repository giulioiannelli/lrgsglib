#!/usr/bin/env python
"""
Complete end-to-end demonstration of entropy analysis for MultiplicativeCascadeGraph.

This script demonstrates the full workflow:
1. Generate eigenvalue data using MCG_SlaplSpect (command-line)
2. Load the eigenvalue arrays
3. Compute entropy observables for each realization
4. Compute ensemble statistics
5. Plot results

Usage:
    # First generate data:
    python MCG_SlaplSpect.py 0.8 0.6 0.6 0.8 0.4 4 \
        --mode eigvals --howmany 0 -p 0.1 --number_of_averages 20 --verbose

    # Then run this analysis:
    python demo_entropy_analysis.py
"""

import pickle as pk
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lrgsglib.utils.lrg.infocomm import compute_entropy_observables_from_eigenvalues
from lrgsglib.graphs.nx import MultiplicativeCascadeGraphNX as MultiplicativeCascadeGraph


def generate_test_data(pflip=0.1, num_realizations=20):
    """
    Generate test eigenvalue data.

    In production, you would run the command-line tool:
    python MCG_SlaplSpect.py 0.8 0.6 0.6 0.8 0.4 5 \
        --mode eigvals --howmany 0 -p 0.1 --number_of_averages 100

    For demo purposes, we generate data programmatically.
    """
    print("=" * 80)
    print("STEP 1: Generate Eigenvalue Data")
    print("=" * 80)
    print(f"\nGenerating {num_realizations} MultiplicativeCascadeGraph realizations...")
    print(f"Parameters: p1=0.8, p2=0.6, p3=0.6, p4=0.8, fraction=0.4, iterations=4")
    print(f"Edge disorder: pflip={pflip}")

    eigvals_list = []

    for i in range(num_realizations):
        # Generate MC graph
        G = MultiplicativeCascadeGraph(
            p1=0.8, p2=0.6, p3=0.6, p4=0.8,
            fraction=0.4,
            iterations=4,  # Small for demo
            variant='exp_clocks',
            pflip=pflip
        )

        # Compute full spectrum (CRITICAL for entropy!)
        G.compute_laplacian_spectrum_weigV()
        eigvals_list.append(G.eigv)

        if i % 5 == 0:
            print(f"  Generated {i+1}/{num_realizations}: N={G.N} nodes")

    print(f"\n✓ Generated {len(eigvals_list)} eigenvalue arrays")
    print(f"  Average system size: {np.mean([len(ev) for ev in eigvals_list]):.1f} nodes")

    return eigvals_list


def compute_entropy_profiles(eigvals_list):
    """
    Compute entropy observables for each realization.
    """
    print("\n" + "=" * 80)
    print("STEP 2: Compute Entropy Observables")
    print("=" * 80)
    print(f"\nComputing entropy profiles for {len(eigvals_list)} realizations...")

    entropy_results = []

    for i, eigvals in enumerate(eigvals_list):
        # Compute entropy observables
        normalized_entropy, specific_heat, variance, time_grid = \
            compute_entropy_observables_from_eigenvalues(
                eigenvalues=eigvals,
                num_nodes=len(eigvals),
                steps=600,
                t1=-2,
                t2=5
            )

        entropy_results.append({
            'entropy': normalized_entropy,
            'specific_heat': specific_heat,
            'variance': variance,
            'time_grid': time_grid
        })

        if i % 5 == 0:
            peak_heat = specific_heat.max()
            print(f"  Processed {i+1}/{len(eigvals_list)}: peak C = {peak_heat:.4f}")

    print(f"\n✓ Computed {len(entropy_results)} entropy profiles")

    return entropy_results


def compute_ensemble_statistics(entropy_results):
    """
    Compute mean and standard deviation over ensemble.
    """
    print("\n" + "=" * 80)
    print("STEP 3: Ensemble Statistics")
    print("=" * 80)

    # Stack results
    entropies = np.array([r['entropy'] for r in entropy_results])
    specific_heats = np.array([r['specific_heat'] for r in entropy_results])
    variances = np.array([r['variance'] for r in entropy_results])
    time_grid = entropy_results[0]['time_grid']

    # Compute statistics
    mean_entropy = entropies.mean(axis=0)
    std_entropy = entropies.std(axis=0)

    mean_specific_heat = specific_heats.mean(axis=0)
    std_specific_heat = specific_heats.std(axis=0)

    mean_variance = variances.mean(axis=0)

    # Find peak in specific heat
    peak_idx = np.argmax(mean_specific_heat)
    peak_time = time_grid[peak_idx]
    peak_value = mean_specific_heat[peak_idx]
    peak_std = std_specific_heat[peak_idx]

    print(f"\nEnsemble averages over {len(entropy_results)} realizations:")
    print(f"  Peak specific heat: C = {peak_value:.4f} ± {peak_std:.4f}")
    print(f"  Peak location: τ = {peak_time:.4f}")
    print(f"  Entropy at peak: S = {mean_entropy[peak_idx]:.4f} ± {std_entropy[peak_idx]:.4f}")
    print(f"  Variance at peak: σ² = {mean_variance[peak_idx]:.4f}")

    return {
        'time_grid': time_grid,
        'mean_entropy': mean_entropy,
        'std_entropy': std_entropy,
        'mean_specific_heat': mean_specific_heat,
        'std_specific_heat': std_specific_heat,
        'mean_variance': mean_variance,
        'peak_time': peak_time,
        'peak_value': peak_value
    }


def plot_results(stats, pflip, save_fig=True):
    """
    Plot entropy profile and specific heat with error bands.
    """
    print("\n" + "=" * 80)
    print("STEP 4: Visualization")
    print("=" * 80)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Plot 1: Entropy profile
    ax = axes[0]
    ax.semilogx(stats['time_grid'], stats['mean_entropy'],
                'b-', linewidth=2, label='Mean')
    ax.fill_between(stats['time_grid'],
                     stats['mean_entropy'] - stats['std_entropy'],
                     stats['mean_entropy'] + stats['std_entropy'],
                     alpha=0.3, label='±1σ')
    ax.axvline(stats['peak_time'], color='r', linestyle='--',
               alpha=0.5, label=f'Peak τ={stats["peak_time"]:.3f}')
    ax.set_xlabel('Time τ', fontsize=14)
    ax.set_ylabel('Normalized Entropy (1 - S)', fontsize=14)
    ax.set_title('Spectral Entropy Profile', fontsize=16, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Plot 2: Specific heat
    ax = axes[1]
    ax.semilogx(stats['time_grid'], stats['mean_specific_heat'],
                'r-', linewidth=2, label='Mean')
    ax.fill_between(stats['time_grid'],
                     stats['mean_specific_heat'] - stats['std_specific_heat'],
                     stats['mean_specific_heat'] + stats['std_specific_heat'],
                     alpha=0.3, label='±1σ')
    ax.axvline(stats['peak_time'], color='r', linestyle='--',
               alpha=0.5, label=f'Peak C={stats["peak_value"]:.3f}')
    ax.set_xlabel('Time τ', fontsize=14)
    ax.set_ylabel('Specific Heat C', fontsize=14)
    ax.set_title('Generalized Specific Heat', fontsize=16, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    # Plot 3: Variance profile
    ax = axes[2]
    ax.semilogx(stats['time_grid'], stats['mean_variance'],
                'g-', linewidth=2, label='Mean variance')
    ax.axvline(stats['peak_time'], color='r', linestyle='--',
               alpha=0.5, label=f'Peak location')
    ax.set_xlabel('Time τ', fontsize=14)
    ax.set_ylabel('Variance σ²', fontsize=14)
    ax.set_title('Spectral Variance', fontsize=16, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    fig.suptitle(f'MultiplicativeCascadeGraph Entropy Analysis (p_flip={pflip})',
                 fontsize=18, fontweight='bold', y=1.02)

    plt.tight_layout()

    if save_fig:
        filename = f'mc_entropy_analysis_pflip_{pflip:.2g}.pdf'
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"\n✓ Figure saved: {filename}")

    plt.show()
    print("✓ Plots generated")


def compare_multiple_pflip_values():
    """
    Compare entropy profiles for different edge disorder values.
    """
    print("\n" + "=" * 80)
    print("BONUS: Comparison Across pflip Values")
    print("=" * 80)
    print("\nComparing spectral entropy for different edge disorder strengths...")

    pflip_values = [0.0, 0.1, 0.2]
    plt.figure(figsize=(14, 6))

    for subplot, (data_type, ylabel, title) in enumerate([
        ('entropy', 'Normalized Entropy (1-S)', 'Entropy vs Edge Disorder'),
        ('specific_heat', 'Specific Heat C', 'Specific Heat vs Edge Disorder')
    ], 1):

        plt.subplot(1, 2, subplot)

        for p in pflip_values:
            print(f"\n  Processing pflip={p}...")

            # Generate data
            eigvals_list = generate_test_data(pflip=p, num_realizations=10)

            # Compute entropy
            results = []
            for eigvals in eigvals_list:
                norm_ent, spec_heat, var, tgrid = \
                    compute_entropy_observables_from_eigenvalues(
                        eigenvalues=eigvals,
                        num_nodes=len(eigvals)
                    )
                if data_type == 'entropy':
                    results.append(norm_ent)
                else:
                    results.append(spec_heat)

            # Average
            mean_curve = np.array(results).mean(axis=0)

            plt.semilogx(tgrid, mean_curve, linewidth=2, label=f'p={p}')

        plt.xlabel('Time τ', fontsize=14)
        plt.ylabel(ylabel, fontsize=14)
        plt.title(title, fontsize=16, fontweight='bold')
        plt.legend(fontsize=12)
        plt.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('mc_entropy_comparison.pdf', dpi=300, bbox_inches='tight')
    print("\n✓ Comparison saved: mc_entropy_comparison.pdf")
    plt.show()


def main():
    """
    Run complete entropy analysis demonstration.
    """
    print("\n" + "=" * 80)
    print("MultiplicativeCascadeGraph Entropy Analysis Demo")
    print("=" * 80)
    print("\nThis demonstrates the complete workflow for computing spectral entropy")
    print("and generalized specific heat from MultiplicativeCascadeGraph eigenvalues.")
    print("=" * 80)

    # Parameters
    pflip = 0.1
    num_realizations = 20

    # Step 1: Generate eigenvalue data
    eigvals_list = generate_test_data(pflip=pflip, num_realizations=num_realizations)

    # Step 2: Compute entropy profiles
    entropy_results = compute_entropy_profiles(eigvals_list)

    # Step 3: Ensemble statistics
    stats = compute_ensemble_statistics(entropy_results)

    # Step 4: Visualization
    plot_results(stats, pflip, save_fig=True)

    # Bonus: Compare different pflip values
    print("\n" + "=" * 80)
    user_input = input("Would you like to see a comparison across pflip values? (y/n): ")
    if user_input.lower() == 'y':
        compare_multiple_pflip_values()

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print("=" * 80)
    print("\nKey findings:")
    print(f"  • Peak specific heat occurs at τ ≈ {stats['peak_time']:.3f}")
    print(f"  • Maximum heat capacity: C ≈ {stats['peak_value']:.3f}")
    print("  • This characterizes the diffusion timescale of the system")
    print("\nFor production runs, use the command-line tool:")
    print("  python MCG_SlaplSpect.py 0.8 0.6 0.6 0.8 0.4 6 \\")
    print("      --mode eigvals --howmany 0 -p 0.1 --number_of_averages 100")
    print("\nThen load and analyze with this script as a template.")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nAnalysis interrupted by user.")
        sys.exit(0)
    except Exception as e:
        print(f"\n\nError during analysis: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
