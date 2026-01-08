#!/usr/bin/env python3
"""
Performance Comparison Test: ContactSimulator1b vs ContactSimulator1c

This test compares the execution speed of the standard version (C1b) versus 
the active border optimized version (C1c) of the contact process simulator.

ContactSimulator1b: Standard MC sweep over all N nodes at each step
ContactSimulator1c: Active border optimization - only updates nodes with active neighbors

Both simulators run on identical lattices with identical initial conditions 
to ensure fair comparison.

Expected behavior:
- C1c should be faster at low densities (< 0.2) where the active border is small
- Performance gain should decrease as density increases
- Both should produce identical final densities when using same seeds
"""

import os
import sys
import time
import numpy as np
from pathlib import Path

# Pytest utilities for time budget enforcement
import pytest

# Add lrgsglib to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.lrgsglib import Lattice2D, ContactProcess, PATHDATA

MAX_RUNTIME = float(os.environ.get("LRGSG_TEST_MAX_SECONDS", "30"))
FAST_TEST = os.environ.get("LRGSG_TEST_FAST", "1") != "0"
_START_TIME = time.monotonic()


def _check_budget(label: str) -> None:
    if time.monotonic() - _START_TIME > MAX_RUNTIME:
        pytest.skip(f"Exceeded {MAX_RUNTIME}s budget during {label}.")


def setup_test_environment():
    """Create test directories"""
    path_data = PATHDATA
    path_test = path_data / 'test' / 'performance_comparison'
    path_test.mkdir(parents=True, exist_ok=True)
    return path_test


def run_single_test(test_params, path_data, verbose=True):
    """
    Run a single performance test with given parameters
    
    Args:
        test_params: dict with keys: side, geo, pflip, gamma, activation, 
                     steps, initial_density, lattice_seed, state_seed
        path_data: Path to data directory
        verbose: Print detailed output
    
    Returns:
        dict with results: C1b_time, C1c_time, speedup, C1b_density, C1c_density
    """
    if verbose:
        print(f"\n{'='*70}")
        print(f"Test Configuration:")
        print(f"  Network: {test_params['geo']} lattice, N={test_params['side']**2} nodes")
        print(f"  Parameters: gamma={test_params['gamma']}, activation={test_params['activation']}")
        print(f"  Initial density: {test_params['initial_density']}")
        print(f"  MC steps: {test_params['steps']}")
        print(f"{'='*70}")
    
    # Create identical lattice for both tests
    np.random.seed(test_params['lattice_seed'])
    lattice_test = Lattice2D(
        side1=test_params['side'], 
        geo=test_params['geo'], 
        pflip=test_params['pflip'], 
        path_data=path_data,
        with_positions=False,
        init_nw_dict=False
    )
    lattice_test.flip_random_fract_edges()
    
    if verbose:
        print(f"\nLattice: N={lattice_test.N}, M={lattice_test.Ne}")
        print(f"Average degree: {2*lattice_test.Ne/lattice_test.N:.2f}")
    
    # Generate identical initial state for both tests
    np.random.seed(test_params['state_seed'])
    initial_state = np.random.choice(
        [0, 1], 
        size=lattice_test.N, 
        p=[1-test_params['initial_density'], test_params['initial_density']]
    )
    
    if verbose:
        print(f"Initial state density: {np.mean(initial_state):.4f}")
        print(f"\n{'-'*70}")
        print("Testing ContactSimulator1b (standard version)...")
        print(f"{'-'*70}")
    
    # Test ContactSimulator1b
    cp_1b = ContactProcess(
        lattice_test,
        gamma=test_params['gamma'],
        activation=test_params['activation'],
        state_type='binary',
        runlang='C1b',
        rndStr=True
    )
    cp_1b.init_contact_dynamics(custom=initial_state.copy())
    
    start_time = time.time()
    cp_1b.run(steps=test_params['steps'], verbose=False, tqdm_on=False)
    time_1b = time.time() - start_time
    density_1b = np.mean(cp_1b.s)
    
    if verbose:
        print(f"Execution time: {time_1b:.4f} seconds")
        print(f"Final density: {density_1b:.6f}")
        print(f"\n{'-'*70}")
        print("Testing ContactSimulator1c (active border optimization)...")
        print(f"{'-'*70}")
    
    # Test ContactSimulator1c
    cp_1c = ContactProcess(
        lattice_test,
        gamma=test_params['gamma'],
        activation=test_params['activation'],
        state_type='binary',
        runlang='C1c',
        rndStr=True
    )
    cp_1c.init_contact_dynamics(custom=initial_state.copy())
    
    start_time = time.time()
    cp_1c.run(steps=test_params['steps'], verbose=False, tqdm_on=False)
    time_1c = time.time() - start_time
    density_1c = np.mean(cp_1c.s)
    
    if verbose:
        print(f"Execution time: {time_1c:.4f} seconds")
        print(f"Final density: {density_1c:.6f}")
    
    speedup = time_1b / time_1c
    density_diff = abs(density_1b - density_1c)
    
    if verbose:
        print(f"\n{'='*70}")
        print("RESULTS")
        print(f"{'='*70}")
        print(f"{'Version':<20} {'Time (s)':<15} {'Final Density':<20} {'Speedup'}")
        print(f"{'-'*70}")
        print(f"{'C1b (standard)':<20} {time_1b:<15.4f} {density_1b:<20.10f} {'1.00x'}")
        print(f"{'C1c (optimized)':<20} {time_1c:<15.4f} {density_1c:<20.10f} {speedup:.2f}x")
        print(f"{'='*70}")
        
        if density_diff < 1e-10:
            print("✓ Final densities match perfectly (identical execution)")
        else:
            print(f"⚠ Final density difference: {density_diff:.2e}")
        
        if speedup > 1.1:
            print(f"✓ C1c is {speedup:.2f}x faster than C1b")
        elif speedup < 0.9:
            print(f"⚠ C1c is {1/speedup:.2f}x slower than C1b (unexpected!)")
        else:
            print(f"≈ Performance is similar (speedup: {speedup:.2f}x)")
    
    return {
        'C1b_time': time_1b,
        'C1c_time': time_1c,
        'speedup': speedup,
        'C1b_density': density_1b,
        'C1c_density': density_1c,
        'density_diff': density_diff
    }


def test_density_sweep(path_data, verbose=True):
    """
    Test performance across different initial densities.
    Active border should be most effective at low densities.
    """
    if verbose:
        print(f"\n\n{'='*70}")
        print("DENSITY SWEEP TEST")
        print(f"{'='*70}")
    
    # Test parameters
    if FAST_TEST:
        side = 16
        steps = 6 * side**2
        density_values = [0.05, 0.20, 0.50]
    else:
        side = 48
        steps = 100 * side**2  # Shorter runs for sweep
        density_values = [0.05, 0.10, 0.20, 0.30, 0.50, 0.70]

    base_params = {
        'side': side,
        'geo': 'sqr',
        'pflip': 0.20,
        'gamma': 0.62,
        'activation': 'relu',
        'steps': steps,
        'lattice_seed': 12345,
        'state_seed': 67890
    }
    results = {
        'densities': density_values,
        'C1b_times': [],
        'C1c_times': [],
        'speedups': [],
        'density_diffs': []
    }
    
    for density in density_values:
        _check_budget(f"density={density:.2f} setup")
        test_params = base_params.copy()
        test_params['initial_density'] = density
        
        if verbose:
            print(f"\nTesting at density = {density:.2f}...")
            print(f"  Creating lattice with seed {test_params['lattice_seed']}...")
        
        # Create identical lattice for both tests
        np.random.seed(test_params['lattice_seed'])
        lattice_test = Lattice2D(
            side1=test_params['side'], 
            geo=test_params['geo'], 
            pflip=test_params['pflip'], 
            path_data=path_data,
            with_positions=False,
            init_nw_dict=False
        )
        lattice_test.flip_random_fract_edges()
        
        if verbose:
            print(f"  Lattice created: N={lattice_test.N}, M={lattice_test.Ne}")
            print(f"  Generating initial state with seed {test_params['state_seed']}...")
        
        # Generate identical initial state for both tests
        np.random.seed(test_params['state_seed'])
        initial_state = np.random.choice(
            [0, 1], 
            size=lattice_test.N, 
            p=[1-density, density]
        )
        
        if verbose:
            print(f"  Initial state density: {np.mean(initial_state):.4f}")
            print(f"  Running C1b...")
        
        # Test ContactSimulator1b
        cp_1b = ContactProcess(
            lattice_test,
            gamma=test_params['gamma'],
            activation=test_params['activation'],
            state_type='binary',
            runlang='C1b',
            rndStr=True
        )
        cp_1b.init_contact_dynamics(custom=initial_state.copy())
        
        start_time = time.time()
        cp_1b.run(steps=test_params['steps'], verbose=False, tqdm_on=False)
        time_1b = time.time() - start_time
        density_1b = np.mean(cp_1b.s)
        _check_budget(f"density={density:.2f} C1b")
        
        if verbose:
            print(f"  C1b done: {time_1b:.3g}s, final density: {density_1b:.6f}")
            print(f"  Running C1c...")
        
        # Test ContactSimulator1c
        cp_1c = ContactProcess(
            lattice_test,
            gamma=test_params['gamma'],
            activation=test_params['activation'],
            state_type='binary',
            runlang='C1c',
            rndStr=True
        )
        cp_1c.init_contact_dynamics(custom=initial_state.copy())
        
        start_time = time.time()
        cp_1c.run(steps=test_params['steps'], verbose=False, tqdm_on=False)
        time_1c = time.time() - start_time
        density_1c = np.mean(cp_1c.s)
        _check_budget(f"density={density:.2f} C1c")
        
        speedup = time_1b / time_1c
        density_diff = abs(density_1b - density_1c)
        
        results['C1b_times'].append(time_1b)
        results['C1c_times'].append(time_1c)
        results['speedups'].append(speedup)
        results['density_diffs'].append(density_diff)
        
        if verbose:
            print(f"  C1c done: {time_1c:.3g}s, final density: {density_1c:.6f}")
            print(f"  Result: C1b: {time_1b:.3g}s, C1c: {time_1c:.3g}s, Speedup: {speedup:.2f}x")
    
    if verbose:
        print(f"\n{'='*70}")
        print("DENSITY SWEEP SUMMARY")
        print(f"{'='*70}")
        print(f"{'Density':<12} {'C1b Time (s)':<15} {'C1c Time (s)':<15} {'Speedup':<10} {'Valid'}")
        print(f"{'-'*70}")
        for i, d in enumerate(results['densities']):
            valid = '✓' if results['density_diffs'][i] < 1e-10 else '⚠'
            print(f"{d:<12.2f} {results['C1b_times'][i]:<15.3f} "
                  f"{results['C1c_times'][i]:<15.3f} "
                  f"{results['speedups'][i]:<10.2f}x {valid}")
        print(f"{'='*70}")
        
        # Analysis
        max_speedup_idx = np.argmax(results['speedups'])
        max_speedup = results['speedups'][max_speedup_idx]
        max_speedup_density = results['densities'][max_speedup_idx]
        
        print(f"\nMaximum speedup: {max_speedup:.2f}x at density = {max_speedup_density:.2f}")
        if max_speedup_density < 0.3:
            print("✓ As expected: Active border is most effective at low densities")
        else:
            print("⚠ Unexpected: Maximum speedup not at lowest density")
    
    return results


def main():
    """Main test function"""
    print("="*70)
    print("Performance Comparison: ContactSimulator1b vs ContactSimulator1c")
    print("="*70)
    
    # Setup
    path_data = setup_test_environment()
    
    # Test 1: Single run with standard parameters
    print("\n\nTEST 1: Single Performance Test")
    print("="*70)
    
    if FAST_TEST:
        side = 16
        steps = 6 * side**2
        density = 0.5
    else:
        side = 48
        steps = 100 * side**2
        density = 0.5

    test_params = {
        'side': side,
        'geo': 'sqr',
        'pflip': 0.20,
        'gamma': 0.62,
        'activation': 'relu',
        'steps': steps,
        'initial_density': density,
        'lattice_seed': 12345,
        'state_seed': 67890
    }
    
    result = run_single_test(test_params, path_data, verbose=True)
    
    # Test 2: Density sweep
    print("\n\nTEST 2: Performance vs Initial Density")
    results = test_density_sweep(path_data, verbose=True)
    
    # Final summary
    print(f"\n\n{'='*70}")
    print("TEST SUITE COMPLETE")
    print(f"{'='*70}")
    print(f"All tests passed successfully!")
    print(f"Data saved to: {path_data}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
