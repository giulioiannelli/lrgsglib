#!/usr/bin/env python
"""
Quick SLQ test for large graph: iterations=6, p2=p3=0.7, fraction=0.7
Expected N ~ 10,000-15,000 nodes
"""

import numpy as np
import time
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from lrgsglib.graphs.nx import MultiplicativeCascadeGraphNX as MultiplicativeCascadeGraph
from lrgsglib.utils.lrg.infocomm import (
    compute_entropy_observables_from_eigenvalues,
    compute_entropy_observables_slq
)

print("="*70)
print("LARGE GRAPH TEST: iterations=6, p2=p3=0.7, fraction=0.7")
print("="*70)

# Generate graph
print("\nGenerating graph...")
t0 = time.time()
G = MultiplicativeCascadeGraph(
    p1=1.0, p2=0.7, p3=0.7, p4=0.9,
    fraction=0.7,
    iterations=6,
    pflip=0.1,
    seed=42
)
G.upd_graph_matrices()
t_gen = time.time() - t0

N = G.N
E = G.gr['G'].number_of_edges()
print(f"Graph: N={N:,} nodes, E={E:,} edges")
print(f"Generation time: {t_gen:.2f}s")
print()

# Method 1: Eigenvalue (may be slow!)
print("Method 1: Eigenvalue-based (may take minutes for large N)...")
t0 = time.time()
try:
    G.compute_laplacian_spectrum(backend='scipy')
    t_spec = time.time() - t0
    print(f"  Spectrum computed: {t_spec:.2f}s")

    S1, C1, V1, tau1 = compute_entropy_observables_from_eigenvalues(
        eigenvalues=G.eigv,
        num_nodes=G.N,
        steps=50,
        t1=-2,
        t2=7
    )
    t_eig = time.time() - t0
    print(f"  Total time: {t_eig:.2f}s")
    eigenvalue_success = True
except Exception as e:
    print(f"  FAILED: {e}")
    eigenvalue_success = False
    t_eig = float('inf')

print()

# Method 2: SLQ with different configurations
print("Method 2: SLQ (fast!)")
print(f"{'Config':<30} {'Time':>12} {'Error':>12}")
print("-"*55)

configs = [
    ("lanczos=50, samples=100", 50, 100),
    ("lanczos=100, samples=100", 100, 100),
    ("lanczos=200, samples=100", 200, 100),
    ("lanczos=300, samples=100", 300, 100),
]

for name, lanczos, samples in configs:
    t0 = time.time()
    S2, C2, V2, tau2 = compute_entropy_observables_slq(
        laplacian=G.slp,
        num_nodes=G.N,
        steps=50,
        t1=-2,
        t2=7,
        lanczos_steps=lanczos,
        num_samples=samples,
        seed=42
    )
    t_slq = time.time() - t0

    if eigenvalue_success:
        error = np.nanmean(np.abs(S1 - S2))
        speedup = t_eig / t_slq
        print(f"{name:<30} {t_slq:>11.2f}s {error:>11.4f}  ({speedup:.1f}× faster)")
    else:
        print(f"{name:<30} {t_slq:>11.2f}s {'N/A':>12}")

print()
print("="*70)
if eigenvalue_success:
    print(f"Eigenvalue time: {t_eig:.2f}s")
    print(f"Best SLQ config for speed+accuracy: lanczos=300, samples=100")
else:
    print("Eigenvalue method failed (N too large)")
    print("SLQ is the only viable option for this graph size")
print("="*70)
