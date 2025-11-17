#!/usr/bin/env python3
"""Quick test to check what mode C1c uses and border statistics"""

import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.lrgsglib import Lattice2D, ContactProcess, PATHDATA

# Small test at very low density
np.random.seed(12345)
lattice = Lattice2D(
    side1=48, 
    geo='sqr', 
    pflip=0.20, 
    path_data=PATHDATA,
    with_positions=False,
    init_nw_dict=False
)
lattice.flip_random_fract_edges()

# 5% initial density
np.random.seed(67890)
initial_state = np.random.choice([0, 1], size=lattice.N, p=[0.95, 0.05])

print(f"Network: N={lattice.N}, M={lattice.Ne}")
print(f"Initial density: {np.mean(initial_state):.4f}")
print(f"\nRunning C1c with verbose output...\n")

cp = ContactProcess(
    lattice,
    gamma=0.62,
    activation='relu',
    state_type='binary',
    runlang='C1c',
    rndStr=True
)
cp.init_contact_dynamics(custom=initial_state.copy())

# Run with stderr visible
cp.run(steps=10 * lattice.N, verbose=True, tqdm_on=False)

print(f"\nFinal density: {np.mean(cp.s):.6f}")
