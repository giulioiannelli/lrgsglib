#!/usr/bin/env python3
"""Test with more steps to see steady-state performance"""

import sys
import time
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.lrgsglib import Lattice2D, ContactProcess, PATHDATA

print("Extended Performance Test\n")

side = 48
steps = 200 * side**2  # Longer run

np.random.seed(12345)
lattice = Lattice2D(side1=side, geo='sqr', pflip=0.20, path_data=PATHDATA,
                   with_positions=False, init_nw_dict=False)
lattice.flip_random_fract_edges()

for initial_density in [0.05, 0.10, 0.20]:
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
    
    # C1c
    cp_1c = ContactProcess(lattice, gamma=0.62, activation='relu',
                          state_type='binary', runlang='C1c', rndStr=True)
    cp_1c.init_contact_dynamics(custom=initial_state.copy())
    t1 = time.time()
    cp_1c.run(steps=steps, verbose=False, tqdm_on=False)
    time_1c = time.time() - t1
    
    speedup = time_1b / time_1c
    print(f"C1b: {time_1b:.3f}s, density: {np.mean(cp_1b.s):.4f}")
    print(f"C1c: {time_1c:.3f}s, density: {np.mean(cp_1c.s):.4f}")
    print(f"Speedup: {speedup:.2f}x {'✓✓' if speedup > 1.2 else '✓' if speedup > 1.05 else '✗'}")
