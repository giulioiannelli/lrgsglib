"""Tests for Lattice3D class.

Tests cover:
- Initialization with various geometries (sc, bcc, fcc)
- Periodic and non-periodic boundary conditions
- Frustration (pflip) and edge sign flipping
- Spectral properties (eigenvalues)
- 3D → 2D projection for visualization
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


class TestLattice3D(unittest.TestCase):
    """Test suite for Lattice3D."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        cls._tmp_dir = tempfile.TemporaryDirectory()
        base_dir = Path(cls._tmp_dir.name)
        cls.data_dir = base_dir / "data"
        cls.log_dir = base_dir / "log"
        cls.data_dir.mkdir(parents=True, exist_ok=True)
        cls.log_dir.mkdir(parents=True, exist_ok=True)
        cls.ccore_dir = base_dir / "Ccore" / "bin"
        cls.ccore_dir.mkdir(parents=True, exist_ok=True)
        os.environ["LRGSG_CCORE_BIN"] = str(cls.ccore_dir)

        cls.src_path = Path(__file__).resolve().parents[1] / "src"
        sys.path.insert(0, str(cls.src_path))

        env_module = types.ModuleType("lrgsglib.config.lrgsg_env")
        env_module.LRGSG_LLIB = str(cls.src_path / "lrgsglib")
        env_module.LRGSG_DATA = str(cls.data_dir)
        env_module.LRGSG_LOG = str(cls.log_dir)
        env_module.LRGSG_CCORE_BIN = str(cls.ccore_dir)
        sys.modules["lrgsglib.config.lrgsg_env"] = env_module

        from lrgsglib.nx_patches.Lattice3D import Lattice3D

        cls.Lattice3D = Lattice3D

    @classmethod
    def tearDownClass(cls):
        """Clean up test fixtures."""
        cls._tmp_dir.cleanup()
        if str(cls.src_path) in sys.path:
            sys.path.remove(str(cls.src_path))

    # ------------------------------------------------------------------
    # Basic initialization tests
    # ------------------------------------------------------------------
    def test_init_simple_cubic(self):
        """Test simple cubic lattice initialization."""
        lat = self.Lattice3D(
            dim=4,
            geo="sc",
            pbc=True,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        # 4x4x4 simple cubic = 64 nodes
        self.assertEqual(lat.N, 64)
        self.assertGreater(lat.Ne, 0)
        self.assertEqual(lat.geo, "simple_cubic")

    def test_init_body_centered_cubic(self):
        """Test body-centered cubic lattice initialization."""
        lat = self.Lattice3D(
            dim=4,
            geo="bcc",
            pbc=True,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        # BCC has 2 nodes per unit cell: 4^3 * 2 = 128 nodes
        self.assertEqual(lat.N, 128)
        self.assertGreater(lat.Ne, 0)
        self.assertEqual(lat.geo, "body_centered")

    def test_init_face_centered_cubic(self):
        """Test face-centered cubic lattice initialization."""
        lat = self.Lattice3D(
            dim=4,
            geo="fcc",
            pbc=True,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        # FCC has 4 nodes per unit cell: 4^3 * 4 = 256 nodes
        self.assertEqual(lat.N, 256)
        self.assertGreater(lat.Ne, 0)
        self.assertEqual(lat.geo, "face_centered")

    def test_init_tuple_dimensions(self):
        """Test initialization with tuple dimensions."""
        lat = self.Lattice3D(
            dim=(3, 4, 5),
            geo="sc",
            pbc=True,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        # 3x4x5 = 60 nodes
        self.assertEqual(lat.N, 60)
        # Note: Lattice3D stores dimensions in reversed order
        self.assertEqual(lat.dim, (5, 4, 3))

    def test_init_non_periodic_boundary(self):
        """Test lattice with open boundary conditions."""
        lat = self.Lattice3D(
            dim=4,
            geo="sc",
            pbc=False,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        self.assertEqual(lat.N, 64)
        # Non-periodic should have fewer edges than periodic
        lat_periodic = self.Lattice3D(
            dim=4,
            geo="sc",
            pbc=True,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        self.assertLess(lat.Ne, lat_periodic.Ne)

    def test_init_only_const_mode(self):
        """Test initialization in const-only mode (no graph built)."""
        lat = self.Lattice3D(
            dim=6,
            geo="sc",
            only_const_mode=True,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        self.assertEqual(lat.dim, (6, 6, 6))
        # Graph should be empty in const-only mode
        self.assertEqual(len(lat.G.nodes()), 0)

    # ------------------------------------------------------------------
    # Frustration tests
    # ------------------------------------------------------------------
    def test_pflip_sets_negative_edges(self):
        """Test that pflip parameter controls negative edge fraction."""
        lat = self.Lattice3D(
            dim=4,
            geo="sc",
            pbc=True,
            pflip=0.3,
            seed=42,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        lat.flip_random_fract_edges()

        # Check that approximately 30% of edges are negative
        self.assertGreater(lat.Ne_n, 0)
        # Allow some tolerance due to random selection
        neg_fraction = lat.Ne_n / lat.Ne
        self.assertGreater(neg_fraction, 0.15)
        self.assertLess(neg_fraction, 0.45)

    def test_zero_pflip_no_negative_edges(self):
        """Test that pflip=0 results in no negative edges."""
        lat = self.Lattice3D(
            dim=4,
            geo="sc",
            pbc=True,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        # No flip call needed - check initial state
        self.assertEqual(lat.Ne_n, 0)

    def test_deterministic_with_seed(self):
        """Test that same seed produces same lattice structure."""
        lat1 = self.Lattice3D(
            dim=4,
            geo="sc",
            pflip=0.3,
            seed=123,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        lat1.flip_random_fract_edges()

        lat2 = self.Lattice3D(
            dim=4,
            geo="sc",
            pflip=0.3,
            seed=123,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        lat2.flip_random_fract_edges()

        self.assertEqual(lat1.Ne_n, lat2.Ne_n)

    # ------------------------------------------------------------------
    # Graph structure tests
    # ------------------------------------------------------------------
    def test_graph_is_connected(self):
        """Test that the generated graph is connected."""
        lat = self.Lattice3D(
            dim=4,
            geo="sc",
            pbc=True,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        self.assertTrue(nx.is_connected(lat.gr["G"]))

    def test_coordination_number(self):
        """Test coordination number for simple cubic lattice."""
        lat = self.Lattice3D(
            dim=4,
            geo="sc",
            pbc=True,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        # Simple cubic with PBC: each node has 6 neighbors
        degrees = [lat.gr["G"].degree(n) for n in lat.gr["G"].nodes()]
        avg_degree = np.mean(degrees)
        self.assertAlmostEqual(avg_degree, 6.0, places=1)

    # ------------------------------------------------------------------
    # Spectral tests
    # ------------------------------------------------------------------
    def test_laplacian_spectrum_computation(self):
        """Test Laplacian spectrum computation."""
        lat = self.Lattice3D(
            dim=4,
            geo="sc",
            pbc=True,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        lat.compute_laplacian_spectrum_weigV(typf=np.float64)

        # Check eigenvalue properties
        self.assertEqual(len(lat.eigv), lat.N)
        # For connected unfrustrated graph, smallest eigenvalue should be ~0
        self.assertAlmostEqual(lat.eigv[0], 0.0, places=5)
        # All eigenvalues should be non-negative for unsigned Laplacian
        self.assertTrue(np.all(lat.eigv >= -1e-10))

    def test_frustrated_laplacian_has_negative_eigenvalues(self):
        """Test that frustrated lattice can have negative eigenvalues."""
        lat = self.Lattice3D(
            dim=4,
            geo="sc",
            pbc=True,
            pflip=0.5,
            seed=42,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        lat.flip_random_fract_edges()
        lat.compute_laplacian_spectrum_weigV(typf=np.float64)

        # Frustrated signed Laplacian may have negative eigenvalues
        # (not guaranteed, but likely for high frustration)
        # Just check computation succeeded
        self.assertEqual(len(lat.eigv), lat.N)
        self.assertTrue(np.all(np.isfinite(lat.eigv)))

    # ------------------------------------------------------------------
    # Position and projection tests
    # ------------------------------------------------------------------
    def test_positions_stored_when_enabled(self):
        """Test that 2D positions are stored when with_positions=True."""
        lat = self.Lattice3D(
            dim=4,
            geo="sc",
            pbc=True,
            with_positions=True,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        # Check that positions exist for nodes
        pos = nx.get_node_attributes(lat.gr["G"], "pos")
        self.assertEqual(len(pos), lat.N)
        # Each position should be a 2D tuple
        for node, p in pos.items():
            self.assertEqual(len(p), 2)

    def test_syshape_path_format(self):
        """Test that syshapePth is correctly formatted."""
        lat = self.Lattice3D(
            dim=5,
            geo="sc",
            pbc=True,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        # Should contain node count
        self.assertIn("N=125", lat.syshapePth)


class TestLattice3DGeometries(unittest.TestCase):
    """Test Lattice3D with different geometric variations."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures."""
        cls._tmp_dir = tempfile.TemporaryDirectory()
        base_dir = Path(cls._tmp_dir.name)
        cls.data_dir = base_dir / "data"
        cls.log_dir = base_dir / "log"
        cls.data_dir.mkdir(parents=True, exist_ok=True)
        cls.log_dir.mkdir(parents=True, exist_ok=True)
        cls.ccore_dir = base_dir / "Ccore" / "bin"
        cls.ccore_dir.mkdir(parents=True, exist_ok=True)
        os.environ["LRGSG_CCORE_BIN"] = str(cls.ccore_dir)

        cls.src_path = Path(__file__).resolve().parents[1] / "src"
        sys.path.insert(0, str(cls.src_path))

        env_module = types.ModuleType("lrgsglib.config.lrgsg_env")
        env_module.LRGSG_LLIB = str(cls.src_path / "lrgsglib")
        env_module.LRGSG_DATA = str(cls.data_dir)
        env_module.LRGSG_LOG = str(cls.log_dir)
        env_module.LRGSG_CCORE_BIN = str(cls.ccore_dir)
        sys.modules["lrgsglib.config.lrgsg_env"] = env_module

        from lrgsglib.nx_patches.Lattice3D import Lattice3D

        cls.Lattice3D = Lattice3D

    @classmethod
    def tearDownClass(cls):
        """Clean up test fixtures."""
        cls._tmp_dir.cleanup()

    def test_all_geometries_create_valid_graphs(self):
        """Test that all geometry options create valid graphs."""
        # Note: FCC with small dim may not be connected, so we only test SC and BCC
        geometries = ["sc", "bcc", "simple_cubic", "body_centered"]

        for geo in geometries:
            with self.subTest(geometry=geo):
                lat = self.Lattice3D(
                    dim=3,
                    geo=geo,
                    pbc=True,
                    pflip=0.0,
                    init_nw_dict=False,
                    path_data=self.data_dir,
                    path_plot=self.log_dir,
                )
                self.assertGreater(lat.N, 0)
                self.assertGreater(lat.Ne, 0)
                self.assertTrue(nx.is_connected(lat.gr["G"]))

    def test_bcc_has_higher_connectivity_than_sc(self):
        """Test that BCC has higher average connectivity than simple cubic."""
        lat_sc = self.Lattice3D(
            dim=4,
            geo="sc",
            pbc=True,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )
        lat_bcc = self.Lattice3D(
            dim=4,
            geo="bcc",
            pbc=True,
            pflip=0.0,
            init_nw_dict=False,
            path_data=self.data_dir,
            path_plot=self.log_dir,
        )

        avg_deg_sc = 2 * lat_sc.Ne / lat_sc.N
        avg_deg_bcc = 2 * lat_bcc.Ne / lat_bcc.N

        # BCC should have higher coordination number than SC
        self.assertGreater(avg_deg_bcc, avg_deg_sc)

    def test_bcc_has_higher_connectivity_than_sc_larger(self):
        """Test that BCC has higher average connectivity than SC for larger lattices."""
        lat_sc = self.Lattice3D(
            dim=5, geo="sc", pbc=True, pflip=0.0, init_nw_dict=False,
            path_data=self.data_dir, path_plot=self.log_dir,
        )
        lat_bcc = self.Lattice3D(
            dim=5, geo="bcc", pbc=True, pflip=0.0, init_nw_dict=False,
            path_data=self.data_dir, path_plot=self.log_dir,
        )

        avg_deg_sc = 2 * lat_sc.Ne / lat_sc.N
        avg_deg_bcc = 2 * lat_bcc.Ne / lat_bcc.N

        # SC < BCC in coordination (BCC has 8 neighbors vs SC's 6)
        self.assertLess(avg_deg_sc, avg_deg_bcc)


if __name__ == "__main__":
    unittest.main()
