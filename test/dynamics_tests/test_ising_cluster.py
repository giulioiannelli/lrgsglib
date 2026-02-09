"""Tests for Ising Wolff and Swendsen-Wang cluster algorithms.

Tests cover all three backends (Python, pybind11, CuPy) for both
cluster algorithms.  Each backend is validated for:
- Basic execution (runs without error, produces valid spin states)
- Output shape (ene, magn, cluster_sizes have correct length)
- Ferromagnetic ordering (low T → aligned state)
- Cluster size tracking (cluster_sizes populated)
- Signed graph support (frustrated lattice produces higher energy)

Cross-backend consistency tests compare Python vs pybind11 (and CuPy
when available) to ensure all implementations agree statistically.
"""

import sys
import tempfile
import types
import unittest
from pathlib import Path

import pytest

np = pytest.importorskip("numpy")
nx = pytest.importorskip("networkx")


# ---------------------------------------------------------------------------
# Check backend availability
# ---------------------------------------------------------------------------

def _has_pybind_cluster():
    """Check if the pybind11 _ising_native module has wolff_sampling."""
    try:
        from lrgsglib.statsys.IsingDynamics.ccore import _ising_native
        return hasattr(_ising_native, "wolff_sampling")
    except ImportError:
        return False


def _has_cupy():
    """Check if CuPy is available."""
    try:
        import cupy  # noqa: F401
        return True
    except ImportError:
        return False


has_pybind = _has_pybind_cluster()
has_cupy = _has_cupy()


# ---------------------------------------------------------------------------
# Fixtures (unittest-style, following existing test patterns)
# ---------------------------------------------------------------------------

class _ClusterTestBase(unittest.TestCase):
    """Shared setup for cluster algorithm tests."""

    @classmethod
    def setUpClass(cls):
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
        cls._tmp_dir.cleanup()

    def _make_lattice(self, side: int = 6, pflip: float = 0.0):
        lat = self.Lattice2D(
            side1=side, side2=side, geo="sqr", pbc=True,
            pflip=pflip, init_nw_dict=False,
            path_data=self.data_dir, path_plot=self.log_dir,
        )
        if pflip > 0:
            lat.flip_random_fract_edges()
        return lat

    def _run_cluster(self, sg, runlang, T=2.0, steps=20, seed=42):
        ising = self.IsingDynamics(
            sg, T=T, steps=steps, runlang=runlang, seed=seed,
        )
        ising.init_ising_dynamics()
        ising.run(tqdm_on=False)
        return ising


# ===================================================================
# Python backend tests
# ===================================================================

class TestPythonWolff(_ClusterTestBase):
    """Tests for the pure-Python Wolff algorithm."""

    def test_runs_without_error(self):
        sg = self._make_lattice(side=6)
        ising = self._run_cluster(sg, "py_wolff", steps=10)
        self.assertTrue(np.all(np.isin(ising.s, (-1, 1))))

    def test_output_shapes(self):
        sg = self._make_lattice(side=6)
        steps = 15
        ising = self._run_cluster(sg, "py_wolff", steps=steps)
        self.assertEqual(len(ising.ene), steps)
        self.assertEqual(len(ising.magn), steps)
        self.assertEqual(len(ising.cluster_sizes), steps)

    def test_energy_is_finite(self):
        sg = self._make_lattice(side=6)
        ising = self._run_cluster(sg, "py_wolff", steps=10)
        for e in ising.ene:
            self.assertTrue(np.isfinite(e))

    def test_magnetization_bounded(self):
        sg = self._make_lattice(side=6)
        ising = self._run_cluster(sg, "py_wolff", steps=10)
        N = sg.N
        for m in ising.magn:
            self.assertGreaterEqual(m, -N)
            self.assertLessEqual(m, N)

    def test_cluster_sizes_positive(self):
        sg = self._make_lattice(side=6)
        ising = self._run_cluster(sg, "py_wolff", steps=10)
        for cs in ising.cluster_sizes:
            self.assertGreater(cs, 0)

    def test_low_T_ferromagnetic_ordering(self):
        sg = self._make_lattice(side=8)
        ising = self._run_cluster(sg, "py_wolff", T=0.5, steps=50, seed=42)
        avg_mag = np.abs(np.mean(ising.s))
        self.assertGreater(avg_mag, 0.8)

    def test_alias_resolves(self):
        """'wolff' alias maps to py_wolff."""
        sg = self._make_lattice(side=4)
        ising = self._run_cluster(sg, "wolff", steps=5)
        self.assertTrue(np.all(np.isin(ising.s, (-1, 1))))

    def test_frustrated_lattice(self):
        sg = self._make_lattice(side=6, pflip=0.3)
        ising = self._run_cluster(sg, "py_wolff", T=1.0, steps=20)
        self.assertTrue(np.all(np.isin(ising.s, (-1, 1))))
        for e in ising.ene:
            self.assertTrue(np.isfinite(e))


class TestPythonSW(_ClusterTestBase):
    """Tests for the pure-Python Swendsen-Wang algorithm."""

    def test_runs_without_error(self):
        sg = self._make_lattice(side=6)
        ising = self._run_cluster(sg, "py_sw", steps=10)
        self.assertTrue(np.all(np.isin(ising.s, (-1, 1))))

    def test_output_shapes(self):
        sg = self._make_lattice(side=6)
        steps = 15
        ising = self._run_cluster(sg, "py_sw", steps=steps)
        self.assertEqual(len(ising.ene), steps)
        self.assertEqual(len(ising.magn), steps)
        self.assertEqual(len(ising.cluster_sizes), steps)

    def test_energy_is_finite(self):
        sg = self._make_lattice(side=6)
        ising = self._run_cluster(sg, "py_sw", steps=10)
        for e in ising.ene:
            self.assertTrue(np.isfinite(e))

    def test_cluster_count_positive(self):
        """SW returns number of clusters (not total flipped)."""
        sg = self._make_lattice(side=6)
        ising = self._run_cluster(sg, "py_sw", steps=10)
        for nc in ising.cluster_sizes:
            self.assertGreater(nc, 0)

    def test_low_T_ferromagnetic_ordering(self):
        sg = self._make_lattice(side=8)
        ising = self._run_cluster(sg, "py_sw", T=0.5, steps=50, seed=42)
        avg_mag = np.abs(np.mean(ising.s))
        self.assertGreater(avg_mag, 0.8)

    def test_alias_resolves(self):
        """'sw' alias maps to py_sw."""
        sg = self._make_lattice(side=4)
        ising = self._run_cluster(sg, "sw", steps=5)
        self.assertTrue(np.all(np.isin(ising.s, (-1, 1))))

    def test_frustrated_lattice(self):
        sg = self._make_lattice(side=6, pflip=0.3)
        ising = self._run_cluster(sg, "py_sw", T=1.0, steps=20)
        self.assertTrue(np.all(np.isin(ising.s, (-1, 1))))


# ===================================================================
# Pybind11 backend tests
# ===================================================================

@pytest.mark.skipif(not has_pybind, reason="pybind11 _ising_native not built")
class TestPybindWolff(_ClusterTestBase):
    """Tests for the pybind11 Wolff algorithm."""

    def test_runs_without_error(self):
        sg = self._make_lattice(side=6)
        ising = self._run_cluster(sg, "pb_wolff", steps=10)
        self.assertTrue(np.all(np.isin(ising.s, (-1, 1))))

    def test_output_shapes(self):
        sg = self._make_lattice(side=6)
        steps = 15
        ising = self._run_cluster(sg, "pb_wolff", steps=steps)
        self.assertEqual(len(ising.ene), steps)
        self.assertEqual(len(ising.magn), steps)
        self.assertEqual(len(ising.cluster_sizes), steps)

    def test_energy_is_finite(self):
        sg = self._make_lattice(side=6)
        ising = self._run_cluster(sg, "pb_wolff", steps=10)
        for e in ising.ene:
            self.assertTrue(np.isfinite(e))

    def test_cluster_sizes_positive(self):
        sg = self._make_lattice(side=6)
        ising = self._run_cluster(sg, "pb_wolff", steps=10)
        for cs in ising.cluster_sizes:
            self.assertGreater(cs, 0)

    def test_low_T_ferromagnetic_ordering(self):
        sg = self._make_lattice(side=8)
        ising = self._run_cluster(sg, "pb_wolff", T=0.5, steps=50, seed=42)
        avg_mag = np.abs(np.mean(ising.s))
        self.assertGreater(avg_mag, 0.8)

    def test_frustrated_lattice(self):
        sg = self._make_lattice(side=6, pflip=0.3)
        ising = self._run_cluster(sg, "pb_wolff", T=1.0, steps=20)
        self.assertTrue(np.all(np.isin(ising.s, (-1, 1))))
        for e in ising.ene:
            self.assertTrue(np.isfinite(e))


@pytest.mark.skipif(not has_pybind, reason="pybind11 _ising_native not built")
class TestPybindSW(_ClusterTestBase):
    """Tests for the pybind11 Swendsen-Wang algorithm."""

    def test_runs_without_error(self):
        sg = self._make_lattice(side=6)
        ising = self._run_cluster(sg, "pb_sw", steps=10)
        self.assertTrue(np.all(np.isin(ising.s, (-1, 1))))

    def test_output_shapes(self):
        sg = self._make_lattice(side=6)
        steps = 15
        ising = self._run_cluster(sg, "pb_sw", steps=steps)
        self.assertEqual(len(ising.ene), steps)
        self.assertEqual(len(ising.magn), steps)
        self.assertEqual(len(ising.cluster_sizes), steps)

    def test_energy_is_finite(self):
        sg = self._make_lattice(side=6)
        ising = self._run_cluster(sg, "pb_sw", steps=10)
        for e in ising.ene:
            self.assertTrue(np.isfinite(e))

    def test_cluster_count_positive(self):
        sg = self._make_lattice(side=6)
        ising = self._run_cluster(sg, "pb_sw", steps=10)
        for nc in ising.cluster_sizes:
            self.assertGreater(nc, 0)

    def test_low_T_ferromagnetic_ordering(self):
        sg = self._make_lattice(side=8)
        ising = self._run_cluster(sg, "pb_sw", T=0.5, steps=50, seed=42)
        avg_mag = np.abs(np.mean(ising.s))
        self.assertGreater(avg_mag, 0.8)

    def test_frustrated_lattice(self):
        sg = self._make_lattice(side=6, pflip=0.3)
        ising = self._run_cluster(sg, "pb_sw", T=1.0, steps=20)
        self.assertTrue(np.all(np.isin(ising.s, (-1, 1))))


# ===================================================================
# CuPy backend tests
# ===================================================================

@pytest.mark.skipif(not has_cupy, reason="CuPy not available")
class TestCuPyWolff(_ClusterTestBase):
    """Tests for the CuPy GPU Wolff algorithm."""

    def test_runs_without_error(self):
        sg = self._make_lattice(side=8)
        ising = self._run_cluster(sg, "cu_wolff", steps=10)
        self.assertTrue(np.all(np.isin(ising.s, (-1, 1))))

    def test_output_shapes(self):
        sg = self._make_lattice(side=8)
        steps = 15
        ising = self._run_cluster(sg, "cu_wolff", steps=steps)
        self.assertEqual(len(ising.ene), steps)
        self.assertEqual(len(ising.magn), steps)
        self.assertEqual(len(ising.cluster_sizes), steps)

    def test_energy_is_finite(self):
        sg = self._make_lattice(side=8)
        ising = self._run_cluster(sg, "cu_wolff", steps=10)
        for e in ising.ene:
            self.assertTrue(np.isfinite(e))

    def test_low_T_ferromagnetic_ordering(self):
        sg = self._make_lattice(side=16)
        ising = self._run_cluster(sg, "cu_wolff", T=0.5, steps=50, seed=42)
        avg_mag = np.abs(np.mean(ising.s))
        self.assertGreater(avg_mag, 0.8)


@pytest.mark.skipif(not has_cupy, reason="CuPy not available")
class TestCuPySW(_ClusterTestBase):
    """Tests for the CuPy GPU Swendsen-Wang algorithm."""

    def test_runs_without_error(self):
        sg = self._make_lattice(side=8)
        ising = self._run_cluster(sg, "cu_sw", steps=10)
        self.assertTrue(np.all(np.isin(ising.s, (-1, 1))))

    def test_output_shapes(self):
        sg = self._make_lattice(side=8)
        steps = 15
        ising = self._run_cluster(sg, "cu_sw", steps=steps)
        self.assertEqual(len(ising.ene), steps)
        self.assertEqual(len(ising.magn), steps)
        self.assertEqual(len(ising.cluster_sizes), steps)

    def test_energy_is_finite(self):
        sg = self._make_lattice(side=8)
        ising = self._run_cluster(sg, "cu_sw", steps=10)
        for e in ising.ene:
            self.assertTrue(np.isfinite(e))

    def test_low_T_ferromagnetic_ordering(self):
        sg = self._make_lattice(side=16)
        ising = self._run_cluster(sg, "cu_sw", T=0.5, steps=50, seed=42)
        avg_mag = np.abs(np.mean(ising.s))
        self.assertGreater(avg_mag, 0.8)


# ===================================================================
# Cross-backend consistency tests
# ===================================================================

class TestCrossBackendConsistency(_ClusterTestBase):
    """Compare outputs across backends for statistical consistency.

    We don't expect exact agreement (different RNGs), but final-state
    magnetization (computed identically from ising.s) should be in
    the same regime.  Energy normalization differs between backends
    (Python: total edge-sum, pybind11: per-site, CuPy: full edge-sum),
    so we compare |<s>| instead.
    """

    def _final_abs_magnetization(self, sg, runlang, T, steps, seed):
        ising = self._run_cluster(sg, runlang, T=T, steps=steps, seed=seed)
        return np.abs(np.mean(ising.s))

    @pytest.mark.skipif(not has_pybind, reason="pybind11 not built")
    def test_wolff_py_vs_pb_low_T(self):
        """Both backends should order at low T."""
        sg = self._make_lattice(side=8)
        T, steps = 0.5, 50

        m_py = self._final_abs_magnetization(sg, "py_wolff", T, steps, seed=42)
        m_pb = self._final_abs_magnetization(sg, "pb_wolff", T, steps, seed=42)

        self.assertGreater(m_py, 0.7, f"py |<s>|={m_py:.3f}")
        self.assertGreater(m_pb, 0.7, f"pb |<s>|={m_pb:.3f}")

    @pytest.mark.skipif(not has_pybind, reason="pybind11 not built")
    def test_sw_py_vs_pb_low_T(self):
        """Both backends should order at low T."""
        sg = self._make_lattice(side=8)
        T, steps = 0.5, 50

        m_py = self._final_abs_magnetization(sg, "py_sw", T, steps, seed=42)
        m_pb = self._final_abs_magnetization(sg, "pb_sw", T, steps, seed=42)

        self.assertGreater(m_py, 0.7, f"py |<s>|={m_py:.3f}")
        self.assertGreater(m_pb, 0.7, f"pb |<s>|={m_pb:.3f}")

    @pytest.mark.skipif(not has_cupy, reason="CuPy not available")
    def test_wolff_py_vs_cu_low_T(self):
        """Both backends should order at low T."""
        sg = self._make_lattice(side=8)
        T, steps = 0.5, 50

        m_py = self._final_abs_magnetization(sg, "py_wolff", T, steps, seed=42)
        m_cu = self._final_abs_magnetization(sg, "cu_wolff", T, steps, seed=42)

        self.assertGreater(m_py, 0.7, f"py |<s>|={m_py:.3f}")
        self.assertGreater(m_cu, 0.7, f"cu |<s>|={m_cu:.3f}")

    @pytest.mark.skipif(not has_cupy, reason="CuPy not available")
    def test_sw_py_vs_cu_low_T(self):
        """Both backends should order at low T."""
        sg = self._make_lattice(side=8)
        T, steps = 0.5, 50

        m_py = self._final_abs_magnetization(sg, "py_sw", T, steps, seed=42)
        m_cu = self._final_abs_magnetization(sg, "cu_sw", T, steps, seed=42)

        self.assertGreater(m_py, 0.7, f"py |<s>|={m_py:.3f}")
        self.assertGreater(m_cu, 0.7, f"cu |<s>|={m_cu:.3f}")


# ===================================================================
# Signed-graph physics tests
# ===================================================================

class TestClusterSignedGraphPhysics(_ClusterTestBase):
    """Verify that cluster algorithms respect signed-graph physics."""

    def test_frustrated_vs_unfrustrated_wolff(self):
        """Frustrated lattice has higher ground-state energy under Wolff."""
        sg_uf = self._make_lattice(side=8, pflip=0.0)
        sg_f = self._make_lattice(side=8, pflip=0.5)

        ising_uf = self._run_cluster(sg_uf, "py_wolff", T=0.5, steps=50, seed=42)
        ising_f = self._run_cluster(sg_f, "py_wolff", T=0.5, steps=50, seed=42)

        e_uf = np.mean(ising_uf.ene[-20:])
        e_f = np.mean(ising_f.ene[-20:])
        self.assertGreater(e_f, e_uf)

    def test_frustrated_vs_unfrustrated_sw(self):
        """Frustrated lattice has higher ground-state energy under SW."""
        sg_uf = self._make_lattice(side=8, pflip=0.0)
        sg_f = self._make_lattice(side=8, pflip=0.5)

        ising_uf = self._run_cluster(sg_uf, "py_sw", T=0.5, steps=50, seed=42)
        ising_f = self._run_cluster(sg_f, "py_sw", T=0.5, steps=50, seed=42)

        e_uf = np.mean(ising_uf.ene[-20:])
        e_f = np.mean(ising_f.ene[-20:])
        self.assertGreater(e_f, e_uf)

    def test_aligned_ferro_stays_aligned_wolff(self):
        """All-up spins on ferro lattice should stay ordered at low T."""
        sg = self._make_lattice(side=6)
        ising = self.IsingDynamics(
            sg, T=0.1, steps=20, runlang="py_wolff",
            ic="custom", seed=42,
        )
        ising.init_ising_dynamics(custom=np.ones(sg.N, dtype=np.int8))
        ising.run(tqdm_on=False)

        # Should remain nearly fully aligned
        avg_mag = np.abs(np.mean(ising.s))
        self.assertGreater(avg_mag, 0.95)

    def test_aligned_ferro_stays_aligned_sw(self):
        """All-up spins on ferro lattice should stay ordered at low T."""
        sg = self._make_lattice(side=6)
        ising = self.IsingDynamics(
            sg, T=0.1, steps=20, runlang="py_sw",
            ic="custom", seed=42,
        )
        ising.init_ising_dynamics(custom=np.ones(sg.N, dtype=np.int8))
        ising.run(tqdm_on=False)

        avg_mag = np.abs(np.mean(ising.s))
        self.assertGreater(avg_mag, 0.95)


# ===================================================================
# Runlang alias resolution tests
# ===================================================================

class TestRunlangAliases(unittest.TestCase):
    """Test runlang alias resolution for cluster algorithms."""

    def test_wolff_alias(self):
        from lrgsglib.statsys.IsingDynamics.IsingDynamics import _resolve_runlang
        self.assertEqual(_resolve_runlang("wolff"), "py_wolff")

    def test_sw_alias(self):
        from lrgsglib.statsys.IsingDynamics.IsingDynamics import _resolve_runlang
        self.assertEqual(_resolve_runlang("sw"), "py_sw")

    def test_canonical_passthrough(self):
        from lrgsglib.statsys.IsingDynamics.IsingDynamics import _resolve_runlang
        self.assertEqual(_resolve_runlang("py_wolff"), "py_wolff")
        self.assertEqual(_resolve_runlang("py_sw"), "py_sw")
        self.assertEqual(_resolve_runlang("pb_wolff"), "pb_wolff")
        self.assertEqual(_resolve_runlang("pb_sw"), "pb_sw")
        self.assertEqual(_resolve_runlang("cu_wolff"), "cu_wolff")
        self.assertEqual(_resolve_runlang("cu_sw"), "cu_sw")


if __name__ == "__main__":
    unittest.main()
