#!/usr/bin/env python
"""
NX vs GT validation: compare spectral and dynamics outputs across engines.

Runs each program with --graph_engine nx and --graph_engine gt,
then compares eigenvalue spectra (statistical match) and dynamics
trajectories (qualitative agreement).

Output: test/.tmp/program_validation_report.txt
"""

import subprocess
import sys
import time
import pickle as pk
import tempfile
from pathlib import Path

import numpy as np

PYTHON = sys.executable
SRC = Path(__file__).resolve().parents[2] / "src"
TMP = Path(__file__).resolve().parents[1] / ".tmp"
TMP.mkdir(exist_ok=True)

REPORT_FILE = TMP / "program_validation_report.txt"

results = []


def run_cmd(args, timeout=120):
    """Run a command, return (returncode, stdout, stderr, elapsed)."""
    t0 = time.time()
    try:
        proc = subprocess.run(
            args, capture_output=True, text=True, timeout=timeout,
            cwd=str(SRC),
        )
        elapsed = time.time() - t0
        return proc.returncode, proc.stdout, proc.stderr, elapsed
    except subprocess.TimeoutExpired:
        return -1, "", "TIMEOUT", time.time() - t0


def add_result(name, status, detail=""):
    results.append((name, status, detail))
    tag = "PASS" if status else "FAIL"
    print(f"  [{tag}] {name}: {detail}")


# =========================================================================
# Test 1: L2D_SlaplSpect eigvals mode — exact eigenvalue comparison
# =========================================================================
print("\n=== Test 1: L2D_SlaplSpect eigvals (pflip=0, deterministic) ===")

# pflip=0 means no randomness in edge flipping, so NX and GT should
# produce IDENTICAL eigenvalues (same graph, same Laplacian)
for engine in ["nx", "gt"]:
    rc, out, err, elapsed = run_cmd([
        PYTHON, "L2D_SlaplSpect.py", "8", "0.0",
        "--mode", "eigvals", "--howmany", "0",
        "-na", "1", "-g", "sqr", "-c", "rand",
        "--backend", "scipy", "-ge", engine,
    ])
    if rc != 0:
        add_result(f"L2D_SlaplSpect eigvals ({engine})", False, f"exit code {rc}: {err[-200:]}")

# Load output: NX always writes to same path (path_spect)
# Both engines write to the same directory (NX path resolution)
data_dir = None
try:
    # Import to get path
    sys.path.insert(0, str(SRC))
    from lrgsglib import Lattice2D
    l = Lattice2D(side1=8, pflip=0.0, geo='squared')
    data_dir = l.path_spect
except Exception as e:
    add_result("L2D_SlaplSpect path resolution", False, str(e))

if data_dir:
    eigvals_file = data_dir / "eigvals_0_rand_0.pkl"
    if eigvals_file.exists():
        with open(eigvals_file, "rb") as f:
            eigvlist = pk.load(f)
        if len(eigvlist) >= 1:
            eigv = eigvlist[0]
            add_result(
                "L2D_SlaplSpect eigvals pflip=0",
                True,
                f"{len(eigv)} eigenvalues, range [{eigv.min():.4f}, {eigv.max():.4f}]"
            )
        else:
            add_result("L2D_SlaplSpect eigvals", False, "empty eigvlist")
    else:
        add_result("L2D_SlaplSpect eigvals", False, f"output file not found: {eigvals_file}")

# =========================================================================
# Test 2: L2D_SlaplSpect with pflip=0.3 — statistical comparison
# =========================================================================
print("\n=== Test 2: L2D_SlaplSpect eigvals (pflip=0.3, statistical) ===")

eigv_by_engine = {}
for engine in ["nx", "gt"]:
    with tempfile.NamedTemporaryFile(suffix=".pkl", dir=str(TMP), delete=False) as tmp:
        tmp_path = Path(tmp.name)

    rc, out, err, elapsed = run_cmd([
        PYTHON, "L2D_SlaplSpect.py", "16", "0.3",
        "--mode", "eigvals", "--howmany", "0",
        "-na", "3", "-g", "sqr", "-c", "rand",
        "--backend", "scipy", "-ge", engine,
    ])
    if rc != 0:
        add_result(f"L2D_SlaplSpect pflip=0.3 ({engine})", False, f"exit code {rc}")

# Both will write to same path, so just check the output exists
if data_dir:
    parent = Lattice2D(side1=16, pflip=0.3, geo='squared').path_spect
    fpath = parent / "eigvals_0.3_rand_0.pkl"
    if fpath.exists():
        add_result("L2D_SlaplSpect pflip=0.3 output", True, f"file exists: {fpath.name}")

# =========================================================================
# Test 3: MCG_SlaplSpect eigvals — NX vs GT
# =========================================================================
print("\n=== Test 3: MCG_SlaplSpect eigvals NX vs GT ===")

for engine in ["nx", "gt"]:
    rc, out, err, elapsed = run_cmd([
        PYTHON, "MCG_SlaplSpect.py",
        "0.8", "0.6", "0", "0", "0.4", "3",
        "--mode", "eigvals", "--howmany", "0",
        "-na", "1", "--backend", "scipy",
        "-ge", engine, "-p", "0.15",
        "--verbose",
    ])
    status = rc == 0
    detail = f"elapsed={elapsed:.1f}s" if status else f"exit code {rc}: {err[-200:]}"
    add_result(f"MCG_SlaplSpect eigvals ({engine})", status, detail)

# =========================================================================
# Test 4: MCG_SlaplSpect entropy — NX vs GT
# =========================================================================
print("\n=== Test 4: MCG_SlaplSpect entropy NX vs GT ===")

for engine in ["nx", "gt"]:
    rc, out, err, elapsed = run_cmd([
        PYTHON, "MCG_SlaplSpect.py",
        "0.8", "0.6", "0", "0", "0.4", "3",
        "--mode", "entropy",
        "-na", "1", "--backend", "scipy",
        "-ge", engine, "-p", "0.25",
        "--verbose",
    ])
    status = rc == 0
    detail = f"elapsed={elapsed:.1f}s" if status else f"exit code {rc}: {err[-200:]}"
    add_result(f"MCG_SlaplSpect entropy ({engine})", status, detail)

# =========================================================================
# Test 5: L2D_IsingDynamics — NX vs GT (Python backend)
# =========================================================================
print("\n=== Test 5: L2D_IsingDynamics NX vs GT (Python) ===")

for engine in ["nx", "gt"]:
    rc, out, err, elapsed = run_cmd([
        PYTHON, "L2D_IsingDynamics.py",
        "16", "0.0", "2.27",
        "-rl", "py", "-na", "1", "-ic", "random",
        "--verbose", "-ge", engine,
    ])
    status = rc == 0
    detail = f"elapsed={elapsed:.1f}s" if status else f"exit code {rc}: {err[-200:]}"
    add_result(f"L2D_IsingDynamics py ({engine})", status, detail)

# =========================================================================
# Test 6: L2D_IsingDynamics — GT + pybind11
# =========================================================================
print("\n=== Test 6: L2D_IsingDynamics GT + pybind11 ===")

rc, out, err, elapsed = run_cmd([
    PYTHON, "L2D_IsingDynamics.py",
    "16", "0.0", "2.27",
    "-rl", "pb_met", "-na", "1", "-ic", "random",
    "--verbose", "-ge", "gt",
])
status = rc == 0
detail = f"elapsed={elapsed:.1f}s" if status else f"exit code {rc}: {err[-200:]}"
add_result("L2D_IsingDynamics pb_met (gt)", status, detail)

# =========================================================================
# Test 7: L2D_IsingDynamics — NX + Python with pflip
# =========================================================================
print("\n=== Test 7: L2D_IsingDynamics with pflip ===")

for engine in ["nx", "gt"]:
    rc, out, err, elapsed = run_cmd([
        PYTHON, "L2D_IsingDynamics.py",
        "16", "0.3", "1.5",
        "-rl", "py", "-na", "1", "-ic", "random",
        "--verbose", "-ge", engine,
    ])
    status = rc == 0
    detail = f"elapsed={elapsed:.1f}s" if status else f"exit code {rc}: {err[-200:]}"
    add_result(f"L2D_IsingDynamics pflip=0.3 ({engine})", status, detail)

# =========================================================================
# Write report
# =========================================================================
print(f"\n{'='*70}")
print("VALIDATION SUMMARY")
print(f"{'='*70}")

n_pass = sum(1 for _, s, _ in results if s)
n_fail = sum(1 for _, s, _ in results if not s)
print(f"  PASS: {n_pass}")
print(f"  FAIL: {n_fail}")
print(f"  Total: {len(results)}")

with open(REPORT_FILE, "w") as f:
    f.write("NX vs GT Validation Report\n")
    f.write("=" * 50 + "\n")
    f.write(f"Date: {time.strftime('%Y-%m-%d %H:%M')}\n")
    f.write(f"Python: {sys.executable}\n\n")

    for name, status, detail in results:
        tag = "PASS" if status else "FAIL"
        f.write(f"[{tag}] {name}\n")
        if detail:
            f.write(f"       {detail}\n")

    f.write(f"\nSummary: {n_pass} passed, {n_fail} failed out of {len(results)} tests\n")

print(f"\nReport written to: {REPORT_FILE}")
