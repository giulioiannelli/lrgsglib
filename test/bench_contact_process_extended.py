#!/usr/bin/env python3
"""Test with more steps to see steady-state performance"""

import os
import sys
import time
import numpy as np
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.lrgsglib import Lattice2D, ContactProcess, PATHDATA

MAX_RUNTIME = float(os.environ.get("LRGSG_TEST_MAX_SECONDS", "30"))
FAST_TEST = os.environ.get("LRGSG_TEST_FAST", "1") != "0"
_START_TIME = time.monotonic()


def _check_budget(label: str) -> None:
    if time.monotonic() - _START_TIME > MAX_RUNTIME:
        pytest.skip(f"Exceeded {MAX_RUNTIME}s budget during {label}.")


print("Extended Performance Test\n")

if FAST_TEST:
    side = 20
    steps = 8 * side**2
    density_values = [0.05]
else:
    side = 48
    steps = 200 * side**2  # Longer run
    density_values = [0.05, 0.10, 0.20]

np.random.seed(12345)
lattice = Lattice2D(side1=side, geo='sqr', pflip=0.20, path_data=PATHDATA,
                   with_positions=False, init_nw_dict=False)
lattice.flip_random_fract_edges()

for initial_density in density_values:
    print(f"\n{'='*60}")
    print(f"Initial density: {initial_density:.2f}")
    print(f"{'='*60}")
    
    np.random.seed(67890)
    initial_state = np.random.choice([0, 1], size=lattice.N, 
                                    p=[1-initial_density, initial_density])
    
    # C1b
    cp_1b = ContactProcess(lattice, gamma=0.62, activation='relu',
                          state_type='binary', runlang='C1b', rndStr=True)
    cp_1b.init_contact_dynamics(custom=initial_state.copy())
    t1 = time.time()
    cp_1b.run(steps=steps, verbose=False, tqdm_on=False)
    time_1b = time.time() - t1
    _check_budget("C1b run")
    
    # C1c
    cp_1c = ContactProcess(lattice, gamma=0.62, activation='relu',
                          state_type='binary', runlang='C1c', rndStr=True)
    cp_1c.init_contact_dynamics(custom=initial_state.copy())
    t1 = time.time()
    cp_1c.run(steps=steps, verbose=False, tqdm_on=False)
    time_1c = time.time() - t1
    _check_budget("C1c run")
    
    speedup = time_1b / time_1c
    print(f"C1b: {time_1b:.3g}s, density: {np.mean(cp_1b.s):.4f}")
    print(f"C1c: {time_1c:.3g}s, density: {np.mean(cp_1c.s):.4f}")
    print(f"Speedup: {speedup:.2f}x {'✓✓' if speedup > 1.2 else '✓' if speedup > 1.05 else '✗'}")
