#!/usr/bin/env python
"""Benchmark MultiplicativeCascadeGraph generation: standard vs exp_clocks."""

import time
import numpy as np
from lrgsglib.nx_patches.MultiplicativeCascade import MultiplicativeCascadeGraph


def benchmark_variant(p1, p2, p3, p4, fraction, iterations, variant, n_runs=10):
    """Benchmark a single variant."""
    times = []

    for i in range(n_runs):
        start = time.time()
        G = MultiplicativeCascadeGraph(
            p1=p1, p2=p2, p3=p3, p4=p4,
            fraction=fraction,
            iterations=iterations,
            variant=variant
        )
        elapsed = time.time() - start
        times.append(elapsed)

        if i == 0:
            n_nodes = G.N
            n_edges = G.Ne

    times = np.array(times)
    return {
        'mean': np.mean(times),
        'std': np.std(times),
        'min': np.min(times),
        'max': np.max(times),
        'n_nodes': n_nodes,
        'n_edges': n_edges
    }


def main():
    """Run benchmarks for different parameter combinations."""
    print("=" * 80)
    print("Benchmarking MultiplicativeCascadeGraph: standard vs exp_clocks")
    print("=" * 80)

    # Test configurations
    configs = [
        {
            'name': 'Small (iterations=4)',
            'p1': 0.8, 'p2': 0.6, 'p3': 0.6, 'p4': 0.8,
            'fraction': 0.4,
            'iterations': 4,
            'n_runs': 20
        },
        {
            'name': 'Medium (iterations=5)',
            'p1': 0.8, 'p2': 0.6, 'p3': 0.6, 'p4': 0.8,
            'fraction': 0.4,
            'iterations': 5,
            'n_runs': 15
        },
        {
            'name': 'Large (iterations=6)',
            'p1': 0.8, 'p2': 0.6, 'p3': 0.6, 'p4': 0.8,
            'fraction': 0.4,
            'iterations': 6,
            'n_runs': 10
        },
        {
            'name': 'Very Large (iterations=7)',
            'p1': 0.8, 'p2': 0.6, 'p3': 0.6, 'p4': 0.8,
            'fraction': 0.4,
            'iterations': 7,
            'n_runs': 5
        },
    ]

    for config in configs:
        print(f"\n{config['name']}")
        print("-" * 80)
        print(f"Parameters: p=({config['p1']}, {config['p2']}, {config['p3']}, {config['p4']}), "
              f"fraction={config['fraction']}, iterations={config['iterations']}")
        print(f"Grid size: {2**(config['iterations']+1)} x {2**(config['iterations']+1)} = "
              f"{(2**(config['iterations']+1))**2} cells")

        # Benchmark standard variant
        print(f"\nTesting 'standard' variant ({config['n_runs']} runs)...")
        standard_result = benchmark_variant(
            config['p1'], config['p2'], config['p3'], config['p4'],
            config['fraction'], config['iterations'],
            variant='standard',
            n_runs=config['n_runs']
        )

        # Benchmark exp_clocks variant
        print(f"Testing 'exp_clocks' variant ({config['n_runs']} runs)...")
        exp_clocks_result = benchmark_variant(
            config['p1'], config['p2'], config['p3'], config['p4'],
            config['fraction'], config['iterations'],
            variant='exp_clocks',
            n_runs=config['n_runs']
        )

        # Print results
        print(f"\nResults:")
        print(f"  Graph size: {standard_result['n_nodes']} nodes, {standard_result['n_edges']} edges")
        print(f"\n  Standard variant:")
        print(f"    Mean time:  {standard_result['mean']*1000:.2f} ± {standard_result['std']*1000:.2f} ms")
        print(f"    Range:      {standard_result['min']*1000:.2f} - {standard_result['max']*1000:.2f} ms")

        print(f"\n  Exp_clocks variant:")
        print(f"    Mean time:  {exp_clocks_result['mean']*1000:.2f} ± {exp_clocks_result['std']*1000:.2f} ms")
        print(f"    Range:      {exp_clocks_result['min']*1000:.2f} - {exp_clocks_result['max']*1000:.2f} ms")

        speedup = standard_result['mean'] / exp_clocks_result['mean']
        print(f"\n  Speedup:    {speedup:.2f}x (exp_clocks is {speedup:.2f}x faster)")

        if speedup > 1:
            print(f"  ✓ exp_clocks is faster by {(speedup-1)*100:.1f}%")
        else:
            print(f"  ✗ standard is faster by {(1/speedup-1)*100:.1f}%")

    print("\n" + "=" * 80)
    print("Benchmark complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
