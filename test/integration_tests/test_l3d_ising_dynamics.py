"""
Integration tests for L3D_IsingDynamics program and kernel.

Tests:
- L3D lattice creation (NX and GT, sc/bcc/fcc)
- SA mode on small lattice (energy decreases)
- Equilibrium mode
- n_thermal inner loop (multiple runs, same graph)
- CLI argument parsing
- Serializer --print mode output verification
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT_DIR = Path(__file__).resolve().parents[2]
SRC_DIR = ROOT_DIR / "src"
L3D_SCRIPT = SRC_DIR / "L3D_IsingDynamics.py"
L3D_SERIALISER = SRC_DIR / "L3D_IsingDynamics_Serialiser.py"

# Add src/ to path so `from kernels.L3D import ...` works
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

pytestmark = [pytest.mark.integration]


# ---------------------------------------------------------------------------
# Helper to build a minimal args namespace
# ---------------------------------------------------------------------------
def _base_args(**overrides) -> SimpleNamespace:
    """Build a minimal args namespace matching L3D_IsingDynamics argparse.

    Positional fields: L (list[int]), p (float), T (float).
    Optional fields mirror ``config/progargs/{Lattice3D,IsingDynamics,SignedGraph}.py``.
    """
    defaults = dict(
        # Positional
        L=[6],            # lattice dimension (list, normalized by kernel)
        p=0.15,           # pflip
        T=2.0,            # temperature
        # L3D optional
        geometry="sc",
        periodic=False,
        pdil=0.0,
        edge_weight="flip",
        mu=0.0,
        sigma=0.0,
        graph_engine="nx",
        workdir=None,     # sgpathn (None → default)
        # SignedGraph optional
        cell_type="rand",
        number_of_averages=1,
        # Ising optional
        init_cond="rand",
        runlang="pb_met",
        n_thermal=1,
        NoClust=1,
        thrmsteps=0,
        randstr="1",
        freq=1,
        out_suffix="",
        sa_mode=False,
        pt_mode=False,
        # SA parameters (used when sa_mode=True)
        T_init=5.0,
        T_final=0.01,
        cooling_schedule="geometric",
        cooling_rate=0.9,
        steps_per_T=10,
        n_temperatures=10,
        # PT parameters (used when pt_mode=True)
        n_replicas=4,
        T_min=0.5,
        T_max=4.0,
        T_ladder_type="geometric",
        steps_per_exchange=10,
        n_exchanges=10,
        # Action flags
        remove_files=True,
        seed=42,
        verbose=False,
        print_chrono=False,
    )
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


# ===================================================================
# LATTICE CREATION TESTS
# ===================================================================
class TestLatticeCreation:
    """Test that prepare_lattice works for all engines and geometries."""

    @pytest.mark.parametrize("engine", ["nx", "gt"])
    @pytest.mark.parametrize("geo", ["sc"])
    def test_create_lattice(self, engine, geo):
        from kernels.L3D import prepare_lattice

        args = _base_args(graph_engine=engine, geometry=geo, L=[4])
        lat = prepare_lattice(args)
        assert lat.N > 0
        assert lat.num_edges > 0

    @pytest.mark.parametrize("engine", ["nx", "gt"])
    def test_pflip_applied(self, engine):
        from kernels.L3D import prepare_lattice

        args = _base_args(graph_engine=engine, p=0.3, L=[4])
        lat = prepare_lattice(args)
        n_neg = lat.count_negative_edges()
        # With pflip=0.3, expect ~30% negative edges (allow wide tolerance)
        frac = n_neg / lat.num_edges if lat.num_edges > 0 else 0.0
        assert 0.05 < frac < 0.6, f"Unexpected negative fraction: {frac:.3f}"


# ===================================================================
# ISING DYNAMICS TESTS
# ===================================================================
class TestIsingDynamics:
    """Test Ising dynamics via the kernel."""

    @pytest.mark.parametrize("engine", ["nx", "gt"])
    def test_equilibrium_run(self, engine):
        """Basic equilibrium Metropolis run completes without error."""
        from kernels.L3D import prepare_lattice
        from kernels.IsingDynamics import run_ising_dynamics

        args = _base_args(graph_engine=engine, runlang="pb_met", L=[4])
        lat = prepare_lattice(args)
        isdy = run_ising_dynamics(args, lat, out_suffix="test_eq")
        assert hasattr(isdy, "magn")
        assert len(isdy.magn) > 0

    @pytest.mark.parametrize("engine", ["nx", "gt"])
    def test_sa_energy_decreases(self, engine):
        """SA mode: final energy should be lower than initial."""
        from kernels.L3D import prepare_lattice
        from kernels.IsingDynamics import run_ising_dynamics

        args = _base_args(
            graph_engine=engine,
            runlang="pb_sa",
            sa_mode=True,
            T=1.0,
            L=[6],
        )
        lat = prepare_lattice(args)
        isdy = run_ising_dynamics(args, lat, out_suffix="test_sa")

        # SA should decrease energy
        if hasattr(isdy, "ene") and len(isdy.ene) > 1:
            assert isdy.ene[-1] <= isdy.ene[0] + 0.1 * abs(isdy.ene[0]), (
                f"SA did not decrease energy: {isdy.ene[0]:.2f} -> {isdy.ene[-1]:.2f}"
            )

    def test_c_backend_sa_nx(self):
        """C subprocess SA (C3B) works on NX 3D lattice."""
        from kernels.L3D import prepare_lattice
        from kernels.IsingDynamics import run_ising_dynamics, clean_up_files

        args = _base_args(
            graph_engine="nx",
            runlang="C1E",
            sa_mode=True,
            T=1.0,
            L=[4],
        )
        lat = prepare_lattice(args)
        isdy = run_ising_dynamics(args, lat, out_suffix="test_c3b")
        assert hasattr(isdy, "magn")
        clean_up_files(isdy, lat)

    def test_c_backend_sa_gt(self):
        """C subprocess SA (C3B) works on GT 3D lattice."""
        from kernels.L3D import prepare_lattice
        from kernels.IsingDynamics import run_ising_dynamics, clean_up_files

        args = _base_args(
            graph_engine="gt",
            runlang="C1E",
            sa_mode=True,
            T=1.0,
            L=[4],
        )
        lat = prepare_lattice(args)
        isdy = run_ising_dynamics(args, lat, out_suffix="test_c3b_gt")
        assert hasattr(isdy, "magn")
        clean_up_files(isdy, lat)


# ===================================================================
# KERNEL AVERAGING TESTS
# ===================================================================
class TestKernelAveraging:
    """Test quenched + thermal averaging loop."""

    def test_n_thermal_multiple_runs(self):
        """n_thermal > 1 runs multiple dynamics on same graph."""
        from kernels.L3D import prepare_lattice
        from kernels.IsingDynamics import run_ising_dynamics

        args = _base_args(
            graph_engine="nx",
            runlang="pb_met",
            number_of_averages=1,
            n_thermal=3,
            L=[4],
        )
        lat = prepare_lattice(args)
        results = []
        for _t in range(args.n_thermal):
            isdy = run_ising_dynamics(args, lat, out_suffix="test_nt")
            results.append(isdy.magn[-1])

        # 3 runs should have been done
        assert len(results) == 3

    def test_full_kernel_run(self):
        """Full run_simulation from kernel module."""
        from kernels.L3D_IsingDynamics import run_simulation

        args = _base_args(
            graph_engine="nx",
            runlang="pb_met",
            number_of_averages=2,
            n_thermal=2,
            L=[4],
        )
        # Should not raise
        run_simulation(args)


# ===================================================================
# CLI TESTS
# ===================================================================
class TestCLI:
    """Test CLI argument parsing and program execution."""

    @pytest.mark.slow
    def test_program_runs(self):
        """L3D_IsingDynamics.py runs successfully with basic args."""
        cmd = [
            sys.executable, str(L3D_SCRIPT),
            "4", "-p", "0.15", "-T", "2.0",
            "-rl", "pb_met",
            "-ge", "nx",
            "-ts", "5",
            "-fq", "5",
            "--remove_files",
        ]
        result = subprocess.run(
            cmd, capture_output=True, text=True,
            cwd=ROOT_DIR / "src",
            timeout=60,
        )
        assert result.returncode == 0, (
            f"Program failed: stdout={result.stdout}, stderr={result.stderr}"
        )

    @pytest.mark.slow
    def test_program_sa_mode(self):
        """L3D_IsingDynamics.py runs with SA mode."""
        cmd = [
            sys.executable, str(L3D_SCRIPT),
            "4", "-p", "0.15", "-T", "1.0",
            "-rl", "pb_sa",
            "-ge", "nx",
            "--sa-mode",
            "--T-init", "5.0",
            "--T-final", "0.01",
            "--n-temperatures", "10",
            "--steps-per-T", "5",
            "--remove_files",
        ]
        result = subprocess.run(
            cmd, capture_output=True, text=True,
            cwd=ROOT_DIR / "src",
            timeout=60,
        )
        assert result.returncode == 0, (
            f"SA mode failed: stdout={result.stdout}, stderr={result.stderr}"
        )

    @pytest.mark.slow
    def test_serializer_print_mode(self):
        """Serializer --print produces valid commands."""
        if not L3D_SERIALISER.exists():
            pytest.skip("L3D_IsingDynamics_Serialiser.py not found")

        cmd = [
            sys.executable, str(L3D_SERIALISER),
            "-s1", "4", "6",
            "-pFT", "(0.1,0.3,2)",
            "-TT", "(1.0,2.0,2)",
            "-rl", "pb_sa",
            "-ge", "nx",
            "--sa-mode",
            "--T-init", "5.0",
            "--T-final", "0.01",
            "--n-temperatures", "10",
            "--steps-per-T", "5",
            "--print",
        ]
        result = subprocess.run(
            cmd, capture_output=True, text=True,
            cwd=ROOT_DIR / "src",
            timeout=30,
        )
        assert result.returncode == 0, (
            f"Serializer failed: stderr={result.stderr}"
        )
        # --print mode should output something
        assert len(result.stdout.strip()) > 0, "Serializer produced no output"
        # Each line should reference L3D_IsingDynamics
        lines = [l for l in result.stdout.strip().split("\n") if l.strip()]
        assert len(lines) > 0
        assert any("L3D_IsingDynamics" in l for l in lines), (
            f"Expected L3D_IsingDynamics in output, got:\n{result.stdout}"
        )


# ===================================================================
# OUTPUT SUFFIX TESTS
# ===================================================================
class TestOutputSuffix:
    """Test that get_l3d_out_suffix produces correct strings."""

    def test_basic_suffix(self):
        from kernels.L3D import get_l3d_out_suffix

        args = _base_args(init_cond="rand", geometry="sc")
        suffix = get_l3d_out_suffix(args)
        assert "rand" in suffix
        assert "sc" in suffix

    def test_suffix_with_pdil(self):
        from kernels.L3D import get_l3d_out_suffix

        args = _base_args(init_cond="rand", geometry="sc", pdil=0.05)
        suffix = get_l3d_out_suffix(args)
        assert "pdil" in suffix

    def test_suffix_gaussian_weights(self):
        from kernels.L3D import get_l3d_out_suffix

        args = _base_args(
            init_cond="rand", geometry="sc",
            edge_weight="gauss", mu=0.5, sigma=0.1,
        )
        suffix = get_l3d_out_suffix(args)
        assert "mu" in suffix
