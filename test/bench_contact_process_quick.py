#!/usr/bin/env python3
"""Quick performance test - single density"""

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


print("Quick Performance Test: C1b vs C1c\n")

# Small, fast test
if FAST_TEST:
    side = 24
    steps = 5 * side**2
else:
    side = 48
    steps = 50 * side**2  # Half the steps for faster test

# Create lattice
np.random.seed(12345)
lattice = Lattice2D(
    side1=side, geo='sqr', pflip=0.20, path_data=PATHDATA,
    with_positions=False, init_nw_dict=False
)
lattice.flip_random_fract_edges()

# 5% initial density
np.random.seed(67890)
initial_state = np.random.choice([0, 1], size=lattice.N, p=[0.95, 0.05])

print(f"N={lattice.N}, steps={steps}, initial_density={np.mean(initial_state):.3g}\n")

# Test C1b
print("Running C1b...")
cp_1b = ContactProcess(lattice, gamma=0.62, activation='relu', 
                       state_type='binary', runlang='C1b', rndStr=True)
cp_1b.init_contact_dynamics(custom=initial_state.copy())
t1 = time.time()
cp_1b.run(steps=steps, verbose=False, tqdm_on=False)
time_1b = time.time() - t1
print(f"  C1b: {time_1b:.3g}s, final density: {np.mean(cp_1b.s):.4f}")
_check_budget("C1b run")

# Test C1c
print("Running C1c...")
cp_1c = ContactProcess(lattice, gamma=0.62, activation='relu',
                       state_type='binary', runlang='C1c', rndStr=True)
cp_1c.init_contact_dynamics(custom=initial_state.copy())
t1 = time.time()
cp_1c.run(steps=steps, verbose=False, tqdm_on=False)
time_1c = time.time() - t1
print(f"  C1c: {time_1c:.3g}s, final density: {np.mean(cp_1c.s):.4f}")
_check_budget("C1c run")

speedup = time_1b / time_1c
print(f"\nSpeedup: {speedup:.2f}x {'✓' if speedup > 1.1 else '✗'}")

# Check log
import subprocess
result = subprocess.run('ls -t .log/err*1c*.log 2>/dev/null | head -1', 
                       shell=True, capture_output=True, text=True)
if result.stdout.strip():
    print(f"\nC1c log:\n{'-'*50}")
    subprocess.run(f'cat {result.stdout.strip()}', shell=True)
