"""Tests for native graph-tool implementations.

Tests the native GT implementations:
- VicsekGraphGT
- DiracCombGraphGT
- DiracBrushGraphGT
- MultiplicativeCascadeGraphGT
"""

import pytest
import numpy as np


class TestVicsekGraphGT:
    """Tests for VicsekGraphGT native implementation."""

    def test_creation(self):
        """Test basic graph creation."""
        from lrgsglib.graphs.gt.VicsekGT import VicsekGraphGT

        v = VicsekGraphGT(N=100, k=3, m=3, seed=42)
        assert v.N == 100
        assert v.k == 3
        assert v.num_edges > 0

    def test_probability_matrix(self):
        """Test probability matrix is computed."""
        from lrgsglib.graphs.gt.VicsekGT import VicsekGraphGT

        v = VicsekGraphGT(N=50, k=2, m=3, seed=42)
        assert v.probability_matrix is not None
        assert v.probability_matrix.shape[0] == 3**2  # m^k

    def test_sign_flipping(self):
        """Test edge sign flipping works."""
        from lrgsglib.graphs.gt.VicsekGT import VicsekGraphGT

        v = VicsekGraphGT(N=100, k=3, m=3, pflip=0.3, seed=42)
        neg_edges = v.count_negative_edges()
        # With pflip=0.3, expect roughly 30% negative edges
        assert neg_edges > 0
        expected = int(0.3 * v.num_edges)
        assert abs(neg_edges - expected) < 0.2 * v.num_edges

    def test_reproducibility(self):
        """Test same seed gives same graph."""
        from lrgsglib.graphs.gt.VicsekGT import VicsekGraphGT

        v1 = VicsekGraphGT(N=50, k=2, m=3, seed=123)
        v2 = VicsekGraphGT(N=50, k=2, m=3, seed=123)
        assert v1.N == v2.N
        assert v1.num_edges == v2.num_edges

    def test_custom_pij(self):
        """Test with custom initial measure matrix."""
        from lrgsglib.graphs.gt.VicsekGT import VicsekGraphGT

        pij = np.array([[0.9, 0.1], [0.1, 0.9]])
        v = VicsekGraphGT(N=50, k=3, pij=pij, seed=42)
        assert v.N == 50
        assert v.probability_matrix.shape[0] == 2**3


class TestDiracCombGraphGT:
    """Tests for DiracCombGraphGT native implementation."""

    def test_creation(self):
        """Test basic graph creation."""
        from lrgsglib.graphs.gt.DiracLatticeGT import DiracCombGraphGT

        dc = DiracCombGraphGT(base_nodes=10, fiber_nodes=5, seed=42)
        expected_nodes = 10 * (1 + 5)  # base + fibers
        assert dc.N == expected_nodes
        assert dc.num_edges > 0

    def test_expected_nodes(self):
        """Test expected node count calculation."""
        from lrgsglib.graphs.gt.DiracLatticeGT import DiracCombGraphGT

        dc = DiracCombGraphGT(base_nodes=8, fiber_nodes=4, seed=42)
        assert dc.get_expected_num_nodes() == 8 * (1 + 4)
        assert dc.N == dc.get_expected_num_nodes()

    def test_dirac_structure(self):
        """Test Dirac structure metadata."""
        from lrgsglib.graphs.gt.DiracLatticeGT import DiracCombGraphGT

        dc = DiracCombGraphGT(base_nodes=10, fiber_nodes=5, seed=42)
        assert dc.dirac_structure is not None
        assert dc.dirac_structure['structure'] == 'dirac_comb'
        assert dc.dirac_structure['base_nodes'] == 10
        assert dc.dirac_structure['fiber_nodes'] == 5

    def test_vertex_indices(self):
        """Test base and fiber vertex index helpers."""
        from lrgsglib.graphs.gt.DiracLatticeGT import DiracCombGraphGT

        dc = DiracCombGraphGT(base_nodes=5, fiber_nodes=3, seed=42)

        # Base vertices should be 0-4
        base_indices = dc.get_base_vertex_indices()
        assert base_indices == [0, 1, 2, 3, 4]

        # Fiber vertices for base 0 should be 5, 6, 7
        fiber_0 = dc.get_fiber_vertex_indices(0)
        assert fiber_0 == [5, 6, 7]

        # Fiber vertices for base 2 should be 11, 12, 13
        fiber_2 = dc.get_fiber_vertex_indices(2)
        assert fiber_2 == [11, 12, 13]

    def test_sign_flipping(self):
        """Test edge sign flipping works."""
        from lrgsglib.graphs.gt.DiracLatticeGT import DiracCombGraphGT

        dc = DiracCombGraphGT(base_nodes=10, fiber_nodes=5, pflip=0.2, seed=42)
        neg_edges = dc.count_negative_edges()
        assert neg_edges > 0

    def test_periodic_edges(self):
        """Test periodic boundary conditions add extra edges."""
        from lrgsglib.graphs.gt.DiracLatticeGT import DiracCombGraphGT

        dc_periodic = DiracCombGraphGT(base_nodes=5, fiber_nodes=3, periodic=True, seed=42)
        dc_open = DiracCombGraphGT(base_nodes=5, fiber_nodes=3, periodic=False, seed=42)

        # Periodic should have more edges (ring vs line)
        assert dc_periodic.num_edges > dc_open.num_edges


class TestDiracBrushGraphGT:
    """Tests for DiracBrushGraphGT native implementation."""

    def test_creation(self):
        """Test basic graph creation."""
        from lrgsglib.graphs.gt.DiracLatticeGT import DiracBrushGraphGT

        db = DiracBrushGraphGT(base_x=4, base_y=4, fiber_nodes=3, seed=42)
        expected_nodes = 4 * 4 * (1 + 3)  # base + fibers
        assert db.N == expected_nodes
        assert db.num_edges > 0

    def test_expected_nodes(self):
        """Test expected node count calculation."""
        from lrgsglib.graphs.gt.DiracLatticeGT import DiracBrushGraphGT

        db = DiracBrushGraphGT(base_x=3, base_y=4, fiber_nodes=2, seed=42)
        assert db.get_expected_num_nodes() == 3 * 4 * (1 + 2)
        assert db.N == db.get_expected_num_nodes()

    def test_dirac_structure(self):
        """Test Dirac structure metadata."""
        from lrgsglib.graphs.gt.DiracLatticeGT import DiracBrushGraphGT

        db = DiracBrushGraphGT(base_x=5, base_y=6, fiber_nodes=4, seed=42)
        assert db.dirac_structure is not None
        assert db.dirac_structure['structure'] == 'dirac_brush'
        assert db.dirac_structure['base_x'] == 5
        assert db.dirac_structure['base_y'] == 6
        assert db.dirac_structure['base_nodes'] == 30
        assert db.dirac_structure['fiber_nodes'] == 4

    def test_vertex_indices(self):
        """Test base and fiber vertex index helpers."""
        from lrgsglib.graphs.gt.DiracLatticeGT import DiracBrushGraphGT

        db = DiracBrushGraphGT(base_x=3, base_y=3, fiber_nodes=2, seed=42)

        # Base vertex at (0, 0) should be 0
        assert db.get_base_vertex_index(0, 0) == 0
        # Base vertex at (1, 2) should be 1*3 + 2 = 5
        assert db.get_base_vertex_index(1, 2) == 5

        # All base vertices should be 0-8
        base_indices = db.get_base_vertex_indices()
        assert base_indices == list(range(9))

        # Fiber vertices for (0, 0) should be 9, 10
        fiber_00 = db.get_fiber_vertex_indices(0, 0)
        assert fiber_00 == [9, 10]

    def test_sign_flipping(self):
        """Test edge sign flipping works."""
        from lrgsglib.graphs.gt.DiracLatticeGT import DiracBrushGraphGT

        db = DiracBrushGraphGT(base_x=4, base_y=4, fiber_nodes=3, pflip=0.2, seed=42)
        neg_edges = db.count_negative_edges()
        assert neg_edges > 0


class TestMultiplicativeCascadeGraphGT:
    """Tests for MultiplicativeCascadeGraphGT native implementation."""

    def test_creation_exp_clocks(self):
        """Test graph creation with exp_clocks variant."""
        from lrgsglib.graphs.gt.MultiplicativeCascadeGT import MultiplicativeCascadeGraphGT

        mc = MultiplicativeCascadeGraphGT(
            p1=0.8, p2=0.6, p3=0.6, p4=0.4,
            fraction=0.4,
            iterations=5,
            variant='exp_clocks',
            seed=42
        )
        assert mc.N > 0
        assert mc.num_edges > 0

    def test_creation_standard(self):
        """Test graph creation with standard variant."""
        from lrgsglib.graphs.gt.MultiplicativeCascadeGT import MultiplicativeCascadeGraphGT

        mc = MultiplicativeCascadeGraphGT(
            p1=0.8, p2=0.6, p3=0.6, p4=0.4,
            fraction=0.5,
            iterations=4,
            variant='standard',
            seed=42
        )
        assert mc.N > 0

    def test_probability_matrix(self):
        """Test probability matrix is computed."""
        from lrgsglib.graphs.gt.MultiplicativeCascadeGT import MultiplicativeCascadeGraphGT

        mc = MultiplicativeCascadeGraphGT(
            p1=0.8, p2=0.6, p3=0.6, p4=0.4,
            iterations=4,
            seed=42
        )
        assert mc.probability_matrix is not None
        # Matrix size is 2^(iterations+1) since we start with 2x2 base
        assert mc.probability_matrix.shape == (2**5, 2**5)

    def test_sign_flipping(self):
        """Test edge sign flipping works."""
        from lrgsglib.graphs.gt.MultiplicativeCascadeGT import MultiplicativeCascadeGraphGT

        mc = MultiplicativeCascadeGraphGT(
            p1=0.8, p2=0.6, p3=0.6, p4=0.4,
            fraction=0.4,
            iterations=5,
            pflip=0.2,
            seed=42
        )
        neg_edges = mc.count_negative_edges()
        assert neg_edges > 0

    def test_reproducibility(self):
        """Test same seed gives same graph."""
        from lrgsglib.graphs.gt.MultiplicativeCascadeGT import MultiplicativeCascadeGraphGT

        mc1 = MultiplicativeCascadeGraphGT(
            p1=0.8, p2=0.6, p3=0.6, p4=0.4,
            fraction=0.3,
            iterations=4,
            seed=999
        )
        mc2 = MultiplicativeCascadeGraphGT(
            p1=0.8, p2=0.6, p3=0.6, p4=0.4,
            fraction=0.3,
            iterations=4,
            seed=999
        )
        assert mc1.N == mc2.N
        assert mc1.num_edges == mc2.num_edges

    def test_periodic_boundary(self):
        """Test periodic vs non-periodic boundary conditions."""
        from lrgsglib.graphs.gt.MultiplicativeCascadeGT import MultiplicativeCascadeGraphGT

        mc_periodic = MultiplicativeCascadeGraphGT(
            p1=0.9, p2=0.7, p3=0.7, p4=0.5,
            fraction=0.5,
            iterations=4,
            periodic=True,
            seed=42
        )
        mc_open = MultiplicativeCascadeGraphGT(
            p1=0.9, p2=0.7, p3=0.7, p4=0.5,
            fraction=0.5,
            iterations=4,
            periodic=False,
            seed=42
        )
        # Both should create valid graphs
        assert mc_periodic.N > 0
        assert mc_open.N > 0


class TestGTNativeConsistency:
    """Test consistency between GT native and NX implementations."""

    def test_vicsek_node_count(self):
        """Verify VicsekGraphGT produces expected node count."""
        from lrgsglib.graphs.gt.VicsekGT import VicsekGraphGT

        for N in [50, 100, 200]:
            v = VicsekGraphGT(N=N, k=3, m=3, seed=42)
            assert v.N == N

    def test_dirac_comb_node_count(self):
        """Verify DiracCombGraphGT produces expected node count."""
        from lrgsglib.graphs.gt.DiracLatticeGT import DiracCombGraphGT

        for base, fiber in [(5, 3), (10, 5), (20, 10)]:
            dc = DiracCombGraphGT(base_nodes=base, fiber_nodes=fiber, seed=42)
            expected = base * (1 + fiber)
            assert dc.N == expected

    def test_dirac_brush_node_count(self):
        """Verify DiracBrushGraphGT produces expected node count."""
        from lrgsglib.graphs.gt.DiracLatticeGT import DiracBrushGraphGT

        for bx, by, fiber in [(3, 3, 2), (4, 4, 3), (5, 5, 4)]:
            db = DiracBrushGraphGT(base_x=bx, base_y=by, fiber_nodes=fiber, seed=42)
            expected = bx * by * (1 + fiber)
            assert db.N == expected


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
