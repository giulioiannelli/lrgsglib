"""Tests for IsingDynamics class.

Tests cover:
- Initialization with various parameters
- Python Metropolis dynamics
- Energy and magnetization calculations
- C backend argument building
- File exports for C backend
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


class TestIsingDynamics(unittest.TestCase):
    """Test suite for IsingDynamics."""

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
        from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph
        from lrgsglib.statsys.IsingDynamics import IsingDynamics

        cls.Lattice2D = Lattice2D
        cls.SignedGraph = SignedGraph
        cls.IsingDynamics = IsingDynamics

    @classmethod
    def tearDownClass(cls):
        """Clean up test fixtures."""
        cls._tmp_dir.cleanup()

    def _make_small_lattice(self, side: int = 4):
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

    def _make_signed_graph(self):
        """Create a small signed graph for testing."""
        G = nx.Graph()
        edges = [
            (0, 1, 1.0),
            (1, 2, -1.0),
            (2, 3, 1.0),
            (3, 0, 1.0),
            (0, 2, 1.0),
            (1, 3, -1.0),
        ]
        for u, v, w in edges:
            G.add_edge(u, v, weight=w)
        return self.SignedGraph(
            G,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )

    # ------------------------------------------------------------------
    # Initialization tests
    # ------------------------------------------------------------------
    def test_init_with_temperature(self):
        """Test initialization with temperature parameter."""
        sg = self._make_small_lattice()
        ising = self.IsingDynamics(sg, T=2.0, steps=10, runlang="py")
        self.assertEqual(ising.T, 2.0)
        self.assertEqual(ising.steps, 10)

    def test_init_with_legacy_eqSTEP(self):
        """Test initialization with legacy eqSTEP parameter."""
        sg = self._make_small_lattice()
        ising = self.IsingDynamics(sg, T=1.5, eqSTEP=20, runlang="py")
        self.assertEqual(ising.steps, 20)
        self.assertEqual(ising.eqSTEP, 20)

    def test_init_with_simref(self):
        """Test initialization with size-normalized time."""
        sg = self._make_small_lattice(side=4)  # 16 nodes
        # Must pass eqSTEP=None to prevent default (20) from overriding simref
        ising = self.IsingDynamics(sg, T=2.0, simref=2.0, eqSTEP=None, runlang="py")
        # simref=2.0 with N=16 should give steps=32
        self.assertEqual(ising.steps, 32)

    def test_init_sa_parameters(self):
        """Test simulated annealing parameter initialization."""
        sg = self._make_small_lattice()
        ising = self.IsingDynamics(
            sg,
            sa_enabled=True,
            T_init=10.0,
            T_final=0.01,
            cooling_schedule="exponential",
            runlang="py",
        )
        self.assertTrue(ising.sa_enabled)
        self.assertEqual(ising.T_init, 10.0)
        self.assertEqual(ising.T_final, 0.01)

    def test_init_pt_parameters(self):
        """Test parallel tempering parameter initialization."""
        sg = self._make_small_lattice()
        ising = self.IsingDynamics(
            sg,
            pt_enabled=True,
            n_replicas=8,
            T_min=0.5,
            T_max=5.0,
            runlang="py",
        )
        self.assertTrue(ising.pt_enabled)
        self.assertEqual(ising.n_replicas, 8)

    # ------------------------------------------------------------------
    # Python dynamics tests
    # ------------------------------------------------------------------
    def test_python_metropolis_runs(self):
        """Test Python Metropolis dynamics runs without error."""
        sg = self._make_small_lattice(side=4)
        ising = self.IsingDynamics(
            sg,
            T=2.0,
            steps=10,
            runlang="py",
            save_magnetization=True,
            seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(tqdm_on=False)

        # Check state is valid (+/-1)
        self.assertTrue(np.all(np.isin(ising.s, (-1, 1))))

    def test_python_dynamics_collects_magnetization(self):
        """Test Python dynamics collects magnetization history."""
        sg = self._make_small_lattice(side=4)  # 16 nodes
        ising = self.IsingDynamics(
            sg,
            T=2.0,
            steps=5,
            runlang="py",
            save_magnetization=True,
            seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(tqdm_on=False)

        # Should have 5 magnetization values
        self.assertEqual(len(ising.magn), 5)
        # Magnetization is stored as SUM of spins (not normalized)
        # For N=16 nodes with spins ±1, the range is [-16, 16]
        N = sg.N
        for m in ising.magn:
            self.assertGreaterEqual(m, -N)
            self.assertLessEqual(m, N)

    def test_python_dynamics_collects_energy(self):
        """Test Python dynamics collects energy history."""
        sg = self._make_small_lattice(side=4)
        ising = self.IsingDynamics(
            sg,
            T=2.0,
            steps=5,
            runlang="py",
            save_magnetization=True,
            seed=42,
        )
        ising.init_ising_dynamics()
        ising.run(tqdm_on=False)

        # Should have 5 energy values
        self.assertEqual(len(ising.ene), 5)
        # Energy should be finite
        for e in ising.ene:
            self.assertTrue(np.isfinite(e))

    def test_low_temperature_ordering(self):
        """Test that low temperature leads to more ordered state."""
        sg = self._make_small_lattice(side=6)

        # Run at low temperature
        ising_low = self.IsingDynamics(
            sg, T=0.5, steps=50, runlang="py", seed=42
        )
        ising_low.init_ising_dynamics()
        ising_low.run(tqdm_on=False)
        mag_low = np.abs(np.mean(ising_low.s))

        # Run at high temperature
        ising_high = self.IsingDynamics(
            sg, T=10.0, steps=50, runlang="py", seed=42
        )
        ising_high.init_ising_dynamics()
        ising_high.run(tqdm_on=False)
        mag_high = np.abs(np.mean(ising_high.s))

        # Low temperature should have higher absolute magnetization
        self.assertGreater(mag_low, mag_high)

    def test_deterministic_with_seed(self):
        """Test that same seed produces same results."""
        sg = self._make_small_lattice(side=4)

        ising1 = self.IsingDynamics(sg, T=2.0, steps=10, runlang="py", seed=123)
        ising1.init_ising_dynamics()
        ising1.run(tqdm_on=False)

        ising2 = self.IsingDynamics(sg, T=2.0, steps=10, runlang="py", seed=123)
        ising2.init_ising_dynamics()
        ising2.run(tqdm_on=False)

        np.testing.assert_array_equal(ising1.s, ising2.s)

    # ------------------------------------------------------------------
    # Energy calculation tests
    # ------------------------------------------------------------------
    def test_calc_full_energy_ferromagnetic(self):
        """Test energy calculation for fully aligned ferromagnetic state."""
        sg = self._make_small_lattice(side=3)
        # Must set ic="custom" to use the custom state
        ising = self.IsingDynamics(sg, T=1.0, ic="custom", runlang="py")
        ising.init_ising_dynamics(custom=np.ones(sg.N, dtype=np.int8))

        energy = ising.calc_full_energy()
        # All spins aligned with positive edges -> minimum energy
        # Energy should be negative for aligned ferromagnet
        self.assertLess(energy, 0)

    def test_calc_full_energy_antiferromagnetic(self):
        """Test energy calculation for antiferromagnetic configuration."""
        sg = self._make_signed_graph()
        ising = self.IsingDynamics(sg, T=1.0, runlang="py")

        # Create alternating pattern
        custom_s = np.array([1, -1, 1, -1], dtype=np.int8)
        ising.init_ising_dynamics(custom=custom_s)

        energy = ising.calc_full_energy()
        # Energy should be finite
        self.assertTrue(np.isfinite(energy))

    # ------------------------------------------------------------------
    # C backend argument building tests
    # ------------------------------------------------------------------
    def test_c_arglist_metropolis(self):
        """Test C argument list for Metropolis variant."""
        sg = self._make_small_lattice(side=4)
        ising = self.IsingDynamics(
            sg, T=2.0, steps=100, runlang="C1b"
        )
        ising.init_ising_dynamics()

        try:
            # Check CbaseName is set correctly
            self.assertEqual(ising.CbaseName, "IsingSimulator1b")

            # Check argument list structure
            arglist = ising._build_c_arglist()
            self.assertIsInstance(arglist, list)
            self.assertGreater(len(arglist), 5)

            # Check temperature is in args
            self.assertIn("2", arglist)
        finally:
            if ising.stderr_fopen and not ising.stderr_fopen.closed:
                ising.stderr_fopen.close()
            ising.remove_run_c_files(remove_stderr=True)
            ising.sg.remove_exported_files()

    def test_c_arglist_sa(self):
        """Test C argument list for simulated annealing variant."""
        sg = self._make_small_lattice(side=4)
        ising = self.IsingDynamics(
            sg,
            sa_enabled=True,
            T_init=10.0,
            T_final=0.01,
            runlang="C3b",
        )
        ising.init_ising_dynamics()

        try:
            self.assertEqual(ising.CbaseName, "IsingSimulator3b")

            arglist = ising._build_c_arglist()
            self.assertIsInstance(arglist, list)

            # Check SA parameters in args
            self.assertIn("10", arglist)  # T_init
            self.assertIn("0.01", arglist)  # T_final
        finally:
            if ising.stderr_fopen and not ising.stderr_fopen.closed:
                ising.stderr_fopen.close()
            ising.remove_run_c_files(remove_stderr=True)
            ising.sg.remove_exported_files()

    def test_c_program_key_validation(self):
        """Test that invalid runlang is rejected."""
        sg = self._make_small_lattice()
        ising = self.IsingDynamics(sg, T=1.0, runlang="C99")

        with self.assertRaises(ValueError):
            ising._c_program_key()

    # ------------------------------------------------------------------
    # File export tests
    # ------------------------------------------------------------------
    def test_init_exports_required_files(self):
        """Test that init_ising_dynamics exports required files for C backend."""
        sg = self._make_small_lattice(side=4)
        ising = self.IsingDynamics(sg, T=2.0, steps=10, runlang="C1b", seed=0)
        ising.init_ising_dynamics()

        try:
            # Check initial state file exists
            self.assertTrue(ising.sfout.exists())

            # Check edge list export
            self.assertTrue(
                hasattr(ising.sg, "path_exp_edgl") and ising.sg.path_exp_edgl.exists()
            )
        finally:
            if ising.stderr_fopen and not ising.stderr_fopen.closed:
                ising.stderr_fopen.close()
            ising.remove_run_c_files(remove_stderr=True)
            ising.sg.remove_exported_files()

    def test_run_c_raises_when_binary_missing(self):
        """Test that running C backend without binary raises FileNotFoundError."""
        sg = self._make_small_lattice(side=4)
        ising = self.IsingDynamics(sg, T=2.0, steps=1, runlang="C1b", seed=0)
        ising.init_ising_dynamics()

        # Point to a non-existent bin directory to simulate missing binary
        from pathlib import Path
        ising._c_bin_dir = Path("/tmp/nonexistent_ccore_bin")
        ising.build_cprogram_command()

        with self.assertRaises(FileNotFoundError):
            ising.run(tqdm_on=False, clean_export=False)

        if ising.stderr_fopen and not ising.stderr_fopen.closed:
            ising.stderr_fopen.close()
        ising.remove_run_c_files(remove_stderr=True)
        ising.sg.remove_exported_files()


class TestIsingDynamicsSignedGraph(unittest.TestCase):
    """Test IsingDynamics with signed (frustrated) graphs."""

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

    def test_frustrated_lattice_dynamics(self):
        """Test dynamics on frustrated lattice (pflip > 0)."""
        lattice = self.Lattice2D(
            side1=6,
            side2=6,
            geo="sqr",
            pbc=True,
            pflip=0.3,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        lattice.flip_random_fract_edges()

        ising = self.IsingDynamics(
            lattice, T=1.0, steps=20, runlang="py", seed=42
        )
        ising.init_ising_dynamics()
        ising.run(tqdm_on=False)

        # Check state is valid
        self.assertTrue(np.all(np.isin(ising.s, (-1, 1))))

    def test_frustrated_vs_unfrustrated_energy(self):
        """Test that frustrated lattice has higher ground state energy."""
        # Unfrustrated lattice
        lattice_uf = self.Lattice2D(
            side1=6, side2=6, geo="sqr", pbc=True, pflip=0.0,
            init_nw_dict=False, path_data=self.data_dir, path_plot=self.log_dir,
        )
        ising_uf = self.IsingDynamics(lattice_uf, T=0.1, steps=50, runlang="py", seed=42)
        ising_uf.init_ising_dynamics()
        ising_uf.run(tqdm_on=False)
        energy_uf = ising_uf.calc_full_energy()

        # Frustrated lattice
        lattice_f = self.Lattice2D(
            side1=6, side2=6, geo="sqr", pbc=True, pflip=0.5,
            init_nw_dict=False, path_data=self.data_dir, path_plot=self.log_dir,
        )
        lattice_f.flip_random_fract_edges()
        ising_f = self.IsingDynamics(lattice_f, T=0.1, steps=50, runlang="py", seed=42)
        ising_f.init_ising_dynamics()
        ising_f.run(tqdm_on=False)
        energy_f = ising_f.calc_full_energy()

        # Frustrated system cannot reach as low energy
        self.assertGreater(energy_f, energy_uf)


if __name__ == "__main__":
    unittest.main()
