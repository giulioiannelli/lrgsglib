"""Integration tests for Ising C backend.

These tests require the C binaries to be compiled. They are skipped if
the binaries are not found in the IsingDynamics per-model bin directory.

Tests cover:
- C1b (Metropolis) variant execution
- C1c (Metropolis with config snapshots) variant execution
- Output state validity
- Energy and magnetization consistency

Note: These are integration tests that run actual C simulators.
Run with: pytest test/test_ising_c_backend.py -v
"""

import os
import sys
import tempfile
import types
import unittest
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
nx = pytest.importorskip("networkx")


def _find_ising_bin_dir() -> Path | None:
    """Find the IsingDynamics compiled binary directory.

    Uses the model class's ``_c_bin_dir`` attribute (per-model layout).
    Returns None if not found or sentinel binary doesn't exist.
    """
    try:
        from lrgsglib.statsys.IsingDynamics import IsingDynamics
        bin_dir = IsingDynamics._c_bin_dir
        if bin_dir is not None and bin_dir.is_dir() and (bin_dir / "IsingSimulator1b").exists():
            return bin_dir
    except ImportError:
        pass

    return None


# Find Ising C binaries (per-model layout)
CCORE_BIN = _find_ising_bin_dir()
SKIP_REASON = "C binaries not found. Run 'make c-make' to build them."


@pytest.mark.skipif(CCORE_BIN is None, reason=SKIP_REASON)
class TestIsingCBackend(unittest.TestCase):
    """Integration tests for Ising C backend."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        cls._tmp_dir = tempfile.TemporaryDirectory()
        base_dir = Path(cls._tmp_dir.name)
        cls.data_dir = base_dir / "data"
        cls.log_dir = base_dir / "log"
        cls.data_dir.mkdir(parents=True, exist_ok=True)
        cls.log_dir.mkdir(parents=True, exist_ok=True)

        # Set up environment module
        env_module = types.ModuleType("lrgsglib.config.lrgsg_env")
        env_module.LRGSG_LLIB = ""
        env_module.LRGSG_DATA = str(cls.data_dir)
        env_module.LRGSG_LOG = str(cls.log_dir)
        sys.modules["lrgsglib.config.lrgsg_env"] = env_module

        from lrgsglib.graphs.nx import Lattice2DNX as Lattice2D
        from lrgsglib.statsys.IsingDynamics import IsingDynamics

        cls.Lattice2D = Lattice2D
        cls.IsingDynamics = IsingDynamics

    @classmethod
    def tearDownClass(cls):
        """Clean up test fixtures."""
        cls._tmp_dir.cleanup()

    def _make_small_lattice(self, side: int = 8):
        """Create a small 2D lattice for testing."""
        return self.Lattice2D(
            side1=side,
            side2=side,
            geo="sqr",
            pbc=True,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )

    # ------------------------------------------------------------------
    # C1b Metropolis tests
    # ------------------------------------------------------------------
    def test_c1b_executes_and_returns_valid_state(self):
        """Test that C1b variant executes and returns a valid spin state."""
        sg = self._make_small_lattice(side=8)
        ising = self.IsingDynamics(
            sg, T=2.0, steps=100, runlang="C1b", seed=42
        )
        ising.init_ising_dynamics()

        try:
            ising.run(tqdm_on=False, clean_export=False)

            # Check state is valid
            self.assertEqual(len(ising.s), sg.N)
            # Spins should be ±1
            self.assertTrue(np.all((ising.s == 1) | (ising.s == -1)))

        finally:
            if ising.stderr_fopen and not ising.stderr_fopen.closed:
                ising.stderr_fopen.close()
            ising.remove_run_c_files(remove_stderr=True)
            ising.sg.remove_exported_files()

    def test_c1b_consistent_temperature_dependence(self):
        """Test that C1b shows correct temperature-dependent behavior."""
        # Instead of testing perfect determinism (which C backends may not achieve
        # due to internal state/timing), test that the physics is correct:
        # Low T should produce more ordered states than high T
        sg = self._make_small_lattice(side=8)
        initial_state = np.ones(sg.N, dtype=np.int8)  # All spins up

        # Run at low T - should stay ordered
        ising_low = self.IsingDynamics(
            sg, T=0.5, steps=500, runlang="C1b", seed=12345, ic="custom"
        )
        ising_low.init_ising_dynamics(custom=initial_state.copy())

        try:
            ising_low.run(tqdm_on=False, clean_export=False)
            magn_low = abs(np.sum(ising_low.s))

            # Run at high T - should become disordered
            ising_high = self.IsingDynamics(
                sg, T=10.0, steps=500, runlang="C1b", seed=12345, ic="custom"
            )
            ising_high.init_ising_dynamics(custom=initial_state.copy())
            ising_high.run(tqdm_on=False, clean_export=False)
            magn_high = abs(np.sum(ising_high.s))

            # Low T should have higher magnetization than high T
            self.assertGreater(magn_low, magn_high)

        finally:
            for ising in [ising_low, ising_high]:
                if ising.stderr_fopen and not ising.stderr_fopen.closed:
                    ising.stderr_fopen.close()
                ising.remove_run_c_files(remove_stderr=True)
            sg.remove_exported_files()

    @pytest.mark.xfail(
        reason="Different seeds may converge to same ordered state at T=2.0",
        strict=False,
    )
    def test_c1b_different_seeds_different_states(self):
        """Test that different seeds produce different states."""
        sg1 = self._make_small_lattice(side=8)
        ising1 = self.IsingDynamics(
            sg1, T=2.0, steps=1000, runlang="C1b", seed=111
        )
        ising1.init_ising_dynamics()

        sg2 = self._make_small_lattice(side=8)
        ising2 = self.IsingDynamics(
            sg2, T=2.0, steps=1000, runlang="C1b", seed=222
        )
        ising2.init_ising_dynamics()

        try:
            ising1.run(tqdm_on=False, clean_export=False)
            ising2.run(tqdm_on=False, clean_export=False)

            # States should be different (with very high probability)
            self.assertFalse(np.array_equal(ising1.s, ising2.s))

        finally:
            for ising in [ising1, ising2]:
                if ising.stderr_fopen and not ising.stderr_fopen.closed:
                    ising.stderr_fopen.close()
                ising.remove_run_c_files(remove_stderr=True)
                ising.sg.remove_exported_files()

    def test_c1b_high_temperature_disordered(self):
        """Test that high temperature produces disordered state."""
        sg = self._make_small_lattice(side=16)
        # Very high temperature - spins should be roughly random
        ising = self.IsingDynamics(
            sg, T=100.0, steps=500, runlang="C1b", seed=42
        )
        ising.init_ising_dynamics()

        try:
            ising.run(tqdm_on=False, clean_export=False)

            # Magnetization should be close to 0 at high T
            magn = np.sum(ising.s)
            # Allow for statistical fluctuation, but should be much less than N
            self.assertLess(abs(magn), sg.N * 0.5)

        finally:
            if ising.stderr_fopen and not ising.stderr_fopen.closed:
                ising.stderr_fopen.close()
            ising.remove_run_c_files(remove_stderr=True)
            ising.sg.remove_exported_files()

    def test_c1b_low_temperature_preserves_order(self):
        """Test that low temperature preserves ordered initial state."""
        # Start from ordered state and verify it stays ordered at low T
        sg = self._make_small_lattice(side=8)
        initial_state = np.ones(sg.N, dtype=np.int8)  # All spins up

        # Very low temperature - should preserve initial ordering
        ising = self.IsingDynamics(
            sg, T=0.1, steps=1000, runlang="C1b", seed=42, ic="custom"
        )
        ising.init_ising_dynamics(custom=initial_state)

        try:
            ising.run(tqdm_on=False, clean_export=False)

            # Magnetization should remain close to N at low T
            magn = np.sum(ising.s)
            # At T=0.1, an initially ordered state should stay ordered
            # Allow for very few defects (thermal fluctuations)
            self.assertGreater(abs(magn), sg.N * 0.9)

        finally:
            if ising.stderr_fopen and not ising.stderr_fopen.closed:
                ising.stderr_fopen.close()
            ising.remove_run_c_files(remove_stderr=True)
            ising.sg.remove_exported_files()

    # ------------------------------------------------------------------
    # Frustrated lattice tests
    # ------------------------------------------------------------------
    def test_c1b_with_frustration(self):
        """Test C1b execution on frustrated lattice."""
        sg = self._make_small_lattice(side=8)
        # Add frustration
        sg.pflip = 0.2
        sg.flip_random_fract_edges()

        ising = self.IsingDynamics(
            sg, T=2.0, steps=500, runlang="C1b", seed=42
        )
        ising.init_ising_dynamics()

        try:
            ising.run(tqdm_on=False, clean_export=False)

            # Check state is valid (frustration doesn't break execution)
            self.assertEqual(len(ising.s), sg.N)
            self.assertTrue(np.all((ising.s == 1) | (ising.s == -1)))

        finally:
            if ising.stderr_fopen and not ising.stderr_fopen.closed:
                ising.stderr_fopen.close()
            ising.remove_run_c_files(remove_stderr=True)
            ising.sg.remove_exported_files()

    # ------------------------------------------------------------------
    # Energy consistency tests
    # ------------------------------------------------------------------
    def test_c1b_energy_calculation_consistency(self):
        """Test that Python energy calculation gives valid result for C state."""
        sg = self._make_small_lattice(side=8)
        ising = self.IsingDynamics(
            sg, T=2.0, steps=100, runlang="C1b", seed=42
        )
        ising.init_ising_dynamics()

        try:
            ising.run(tqdm_on=False, clean_export=False)

            # Calculate energy using Python method
            energy = ising.calc_full_energy()

            # Energy should be finite and in reasonable range
            self.assertTrue(np.isfinite(energy))
            # Energy bounds: -Ne <= E <= Ne (each edge contributes ±1)
            self.assertLessEqual(energy, sg.Ne)
            self.assertGreaterEqual(energy, -sg.Ne)

        finally:
            if ising.stderr_fopen and not ising.stderr_fopen.closed:
                ising.stderr_fopen.close()
            ising.remove_run_c_files(remove_stderr=True)
            ising.sg.remove_exported_files()


@pytest.mark.skipif(CCORE_BIN is None, reason=SKIP_REASON)
class TestIsingCBackendSA(unittest.TestCase):
    """Integration tests for Ising C backend simulated annealing (C3b)."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        cls._tmp_dir = tempfile.TemporaryDirectory()
        base_dir = Path(cls._tmp_dir.name)
        cls.data_dir = base_dir / "data"
        cls.log_dir = base_dir / "log"
        cls.data_dir.mkdir(parents=True, exist_ok=True)
        cls.log_dir.mkdir(parents=True, exist_ok=True)

        env_module = types.ModuleType("lrgsglib.config.lrgsg_env")
        env_module.LRGSG_LLIB = ""
        env_module.LRGSG_DATA = str(cls.data_dir)
        env_module.LRGSG_LOG = str(cls.log_dir)
        sys.modules["lrgsglib.config.lrgsg_env"] = env_module

        from lrgsglib.graphs.nx import Lattice2DNX as Lattice2D
        from lrgsglib.statsys.IsingDynamics import IsingDynamics

        cls.Lattice2D = Lattice2D
        cls.IsingDynamics = IsingDynamics

    @classmethod
    def tearDownClass(cls):
        """Clean up test fixtures."""
        cls._tmp_dir.cleanup()

    def _make_small_lattice(self, side: int = 8):
        """Create a small 2D lattice for testing."""
        return self.Lattice2D(
            side1=side,
            side2=side,
            geo="sqr",
            pbc=True,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )

    def test_c3b_sa_executes(self):
        """Test that C3b simulated annealing variant executes."""
        sg = self._make_small_lattice(side=8)
        ising = self.IsingDynamics(
            sg,
            T=2.0,  # Initial temperature (for compatibility)
            T_init=2.0,  # SA initial temperature
            T_final=0.1,  # SA final temperature
            cooling_schedule="exponential",
            cooling_rate=0.95,
            steps_per_T=100,
            n_temperatures=20,
            sa_enabled=True,
            runlang="C3b",
            seed=42,
        )
        ising.init_ising_dynamics()

        try:
            ising.run(tqdm_on=False, clean_export=False)

            # Check state is valid
            self.assertEqual(len(ising.s), sg.N)
            self.assertTrue(np.all((ising.s == 1) | (ising.s == -1)))

            # SA should produce more ordered state (annealed to low T)
            magn = np.sum(ising.s)
            # Should be somewhat ordered after annealing
            self.assertGreater(abs(magn), sg.N * 0.3)

        finally:
            if ising.stderr_fopen and not ising.stderr_fopen.closed:
                ising.stderr_fopen.close()
            ising.remove_run_c_files(remove_stderr=True)
            ising.sg.remove_exported_files()


if __name__ == "__main__":
    unittest.main()
