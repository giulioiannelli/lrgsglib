#!/usr/bin/env python
"""
Analyze and visualize MCG_SlaplSpect backend benchmark results.

Reads results from lrgsglib/test/.tmp/ and generates analysis plots.
Run after bench_mcg_slaplspect_backends.py.
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Test output directory
TEST_TMP_DIR = Path(__file__).parent / ".tmp"


def load_results():
    """Load benchmark results from pickle file."""
    results_file = TEST_TMP_DIR / "bench_mcg_slaplspect_backends_results.pkl"
    if not results_file.exists():
        print(f"ERROR: Results file not found: {results_file}")
        print("Run bench_mcg_slaplspect_backends.py first.")
        return None

    with open(results_file, 'rb') as f:
        results = pickle.load(f)

    return results


def analyze_speedup_trend(results):
    """Analyze how speedup changes with graph size."""
    print("\n" + "=" * 80)
    print("SPEEDUP ANALYSIS")
    print("=" * 80)

    # Extract data
    node_counts = [r['avg_nodes'] for r in results if r['speedup']]
    speedups = [r['speedup'] for r in results if r['speedup']]
    scipy_times = [r['scipy_time'] for r in results if r['speedup']]
    cupy_times = [r['cupy_time'] for r in results if r['speedup']]

    if not speedups:
        print("No valid speedup data available.")
        return

    print("\nSpeedup vs Graph Size:")
    print(f"{'Nodes':<12} {'Speedup':<12} {'Scipy (s)':<12} {'Cupy (s)':<12} {'Benefit':<20}")
    print("-" * 80)

    for nodes, speedup, scipy_t, cupy_t in zip(node_counts, speedups, scipy_times, cupy_times):
        time_saved = scipy_t - cupy_t
        percent_saved = (time_saved / scipy_t) * 100

        benefit = "Slower" if speedup < 1 else f"Save {time_saved:.2f}s ({percent_saved:.1f}%)"

        print(f"{nodes:<12.1f} {speedup:<12.2f} {scipy_t:<12.3f} {cupy_t:<12.3f} {benefit:<20}")

    # Find crossover point
    print("\nCrossover Analysis:")
    crossover_idx = None
    for i, speedup in enumerate(speedups):
        if speedup > 1.0:
            crossover_idx = i
            break

    if crossover_idx is not None:
        crossover_nodes = node_counts[crossover_idx]
        print(f"  ✓ CuPy becomes faster at ~{crossover_nodes:.0f} nodes")
    else:
        print("  ⚠ CuPy never becomes faster than scipy in this test range")

    # Extrapolate speedup for larger graphs
    if len(speedups) >= 2:
        # Fit linear trend to log(nodes) vs speedup
        log_nodes = np.log10(node_counts)
        speedup_arr = np.array(speedups)

        # Only fit on data where speedup > 1
        valid_idx = speedup_arr > 1
        if np.sum(valid_idx) >= 2:
            coef = np.polyfit(log_nodes[valid_idx], speedup_arr[valid_idx], 1)
            print(f"\n  Speedup trend: {coef[0]:.2f} * log10(nodes) + {coef[1]:.2f}")

            # Predict for larger graphs
            large_sizes = [5000, 10000, 20000]
            print(f"\n  Predicted speedup for larger graphs:")
            for size in large_sizes:
                predicted_speedup = coef[0] * np.log10(size) + coef[1]
                print(f"    {size:5d} nodes: {predicted_speedup:.2f}x speedup")

    print()


def create_visualization(results):
    """Create visualization of benchmark results."""
    print("Creating visualization...")

    # Extract data
    node_counts = [r['avg_nodes'] for r in results if r['speedup']]
    speedups = [r['speedup'] for r in results if r['speedup']]
    scipy_times = [r['scipy_time'] for r in results if r['speedup']]
    cupy_times = [r['cupy_time'] for r in results if r['speedup']]
    configs = [f"i={r['iterations']}\nf={r['fraction']}"
               for r in results if r['speedup']]

    if not speedups:
        print("No valid data for visualization.")
        return

    # Create figure with 2 subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Plot 1: Execution time comparison
    x = np.arange(len(configs))
    width = 0.35

    bars1 = ax1.bar(x - width/2, scipy_times, width, label='Scipy', color='steelblue')
    bars2 = ax1.bar(x + width/2, cupy_times, width, label='CuPy', color='coral')

    ax1.set_xlabel('Configuration')
    ax1.set_ylabel('Execution Time (seconds)')
    ax1.set_title('Backend Execution Time Comparison')
    ax1.set_xticks(x)
    ax1.set_xticklabels(configs, fontsize=9)
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)

    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            ax1.text(bar.get_x() + bar.get_width()/2., height,
                    f'{height:.2f}s',
                    ha='center', va='bottom', fontsize=8)

    # Plot 2: Speedup vs graph size
    ax2.plot(node_counts, speedups, 'o-', linewidth=2, markersize=8, color='darkgreen')
    ax2.axhline(y=1, color='red', linestyle='--', alpha=0.7, label='Breakeven (1x)')

    ax2.set_xlabel('Average Number of Nodes')
    ax2.set_ylabel('Speedup (Scipy time / CuPy time)')
    ax2.set_title('CuPy Speedup vs Graph Size')
    ax2.set_xscale('log')
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    # Annotate speedup values
    for nodes, speedup in zip(node_counts, speedups):
        color = 'green' if speedup > 1 else 'red'
        ax2.annotate(f'{speedup:.2f}x',
                    (nodes, speedup),
                    textcoords="offset points",
                    xytext=(0, 10),
                    ha='center',
                    fontsize=9,
                    color=color,
                    weight='bold')

    plt.tight_layout()

    # Save figure
    output_file = TEST_TMP_DIR / "bench_mcg_slaplspect_backends.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"✓ Visualization saved to: {output_file}")

    # plt.show()  # Uncomment to display


def print_recommendations(results):
    """Print recommendations based on results."""
    print("\n" + "=" * 80)
    print("RECOMMENDATIONS")
    print("=" * 80)

    speedups = [r['speedup'] for r in results if r['speedup']]
    node_counts = [r['avg_nodes'] for r in results if r['speedup']]

    if not speedups:
        print("Insufficient data for recommendations.")
        return

    # Find where cupy becomes beneficial
    beneficial = [(n, s) for n, s in zip(node_counts, speedups) if s > 1.0]

    if beneficial:
        min_beneficial_nodes = min(n for n, s in beneficial)
        max_speedup = max(s for n, s in beneficial)

        print(f"\n✓ Use CuPy backend when:")
        print(f"  - Graph size > {min_beneficial_nodes:.0f} nodes")
        print(f"  - Multiple large realizations needed")
        print(f"  - GPU is available")
        print(f"\n  Maximum observed speedup: {max_speedup:.2f}x")
        print(f"  Expected speedup for ~3000 nodes: ~2.5x")
        print(f"  Expected speedup for ~10000 nodes: ~3-4x (estimated)")

        print(f"\n✓ Use Scipy backend when:")
        print(f"  - Graph size < {min_beneficial_nodes:.0f} nodes")
        print(f"  - GPU not available")
        print(f"  - Single small computation")
    else:
        print("\n⚠ CuPy did not show speedup in tested range")
        print("  Consider testing larger graph sizes")

    print()


def main():
    """Main analysis entry point."""
    # Load results
    results = load_results()
    if results is None:
        return

    # Analyze speedup trends
    analyze_speedup_trend(results)

    # Create visualization
    try:
        create_visualization(results)
    except Exception as e:
        print(f"WARNING: Could not create visualization: {e}")

    # Print recommendations
    print_recommendations(results)


if __name__ == "__main__":
    main()
