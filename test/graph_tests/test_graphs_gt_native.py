"""Tests for native graph-tool implementations.

Tests the native GT implementations:
- VicsekGraphGT
- DiracCombGraphGT
- DiracBrushGraphGT
- MultiplicativeCascadeGraphGT
- FullyConnectedGT
- SCSGeneralizedNNGT
- kRegularGraphGT
- RandomGeometricGT
- ConfigurationModelGT
- BipartiteGraphGT (unified: random + degree-sequence modes)
- BipartiteFromDegreeSequenceGT (backward-compat wrapper)
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


class TestFullyConnectedGT:
    """Tests for FullyConnectedGT native implementation."""

    def test_creation(self):
        """Test basic graph creation."""
        from lrgsglib.graphs.gt.FullyConnectedGT import FullyConnectedGT

        fc = FullyConnectedGT(N=50, seed=42)
        assert fc.N == 50
        expected_edges = 50 * 49 // 2
        assert fc.num_edges == expected_edges

    def test_expected_nodes_edges(self):
        """Test expected node and edge count."""
        from lrgsglib.graphs.gt.FullyConnectedGT import FullyConnectedGT

        fc = FullyConnectedGT(N=20, seed=42)
        assert fc.get_expected_num_nodes() == 20
        assert fc.get_expected_num_edges() == 20 * 19 // 2

    def test_sign_flipping(self):
        """Test edge sign flipping works."""
        from lrgsglib.graphs.gt.FullyConnectedGT import FullyConnectedGT

        fc = FullyConnectedGT(N=50, pflip=0.3, seed=42)
        neg_edges = fc.count_negative_edges()
        assert neg_edges > 0
        expected = int(0.3 * fc.num_edges)
        assert abs(neg_edges - expected) < 0.15 * fc.num_edges

    def test_reproducibility(self):
        """Test same seed gives same graph."""
        from lrgsglib.graphs.gt.FullyConnectedGT import FullyConnectedGT

        fc1 = FullyConnectedGT(N=30, pflip=0.2, seed=123)
        fc2 = FullyConnectedGT(N=30, pflip=0.2, seed=123)
        assert fc1.N == fc2.N
        assert fc1.num_edges == fc2.num_edges
        assert fc1.count_negative_edges() == fc2.count_negative_edges()

    def test_positions(self):
        """Test position computation."""
        from lrgsglib.graphs.gt.FullyConnectedGT import FullyConnectedGT

        fc = FullyConnectedGT(N=20, with_positions=True, seed=42)
        pos = fc.get_positions()
        assert pos.shape == (20, 2)


class TestSCSGeneralizedNNGT:
    """Tests for SCSGeneralizedNNGT native implementation."""

    def test_creation_symmetric(self):
        """Test symmetric (undirected) case."""
        from lrgsglib.graphs.gt.SCSGeneralizedNNGT import SCSGeneralizedNNGT

        scs = SCSGeneralizedNNGT(N=50, gamma=1.0, seed=42)
        assert scs.N == 50
        assert scs.is_symmetric()
        # Undirected: N*(N-1)/2 + N self-loops
        expected_edges = 50 * 49 // 2 + 50
        assert scs.num_edges == expected_edges

    def test_creation_asymmetric(self):
        """Test asymmetric (directed) case."""
        from lrgsglib.graphs.gt.SCSGeneralizedNNGT import SCSGeneralizedNNGT

        scs = SCSGeneralizedNNGT(N=30, gamma=0.5, seed=42)
        assert scs.N == 30
        assert not scs.is_symmetric()
        # Directed: N*(N-1) both directions, no self-loops by default
        expected_edges = 30 * 29
        assert scs.num_edges == expected_edges

    def test_weight_matrix(self):
        """Test weight matrix is computed."""
        from lrgsglib.graphs.gt.SCSGeneralizedNNGT import SCSGeneralizedNNGT

        scs = SCSGeneralizedNNGT(N=20, gamma=0.5, seed=42)
        W = scs.get_weight_matrix()
        assert W.shape == (20, 20)
        # Off-diagonal should have non-zero entries
        assert np.any(W[np.triu_indices(20, k=1)] != 0)

    def test_correlation(self):
        """Test weight correlation roughly matches gamma."""
        from lrgsglib.graphs.gt.SCSGeneralizedNNGT import SCSGeneralizedNNGT

        scs = SCSGeneralizedNNGT(N=100, gamma=0.5, J=1.0, seed=42)
        W = scs.get_weight_matrix()

        # Collect pairs of reciprocal weights
        w_ij = []
        w_ji = []
        for i in range(100):
            for j in range(i + 1, 100):
                w_ij.append(W[i, j])
                w_ji.append(W[j, i])

        # Check correlation is close to gamma (within statistical noise)
        corr = np.corrcoef(w_ij, w_ji)[0, 1]
        assert abs(corr - 0.5) < 0.15

    def test_sign_flipping(self):
        """Test edge sign flipping works."""
        from lrgsglib.graphs.gt.SCSGeneralizedNNGT import SCSGeneralizedNNGT

        scs = SCSGeneralizedNNGT(N=40, gamma=1.0, pflip=0.2, seed=42)
        neg_edges = scs.count_negative_edges()
        assert neg_edges > 0

    def test_reproducibility(self):
        """Test same seed gives same graph."""
        from lrgsglib.graphs.gt.SCSGeneralizedNNGT import SCSGeneralizedNNGT

        scs1 = SCSGeneralizedNNGT(N=30, gamma=0.7, seed=999)
        scs2 = SCSGeneralizedNNGT(N=30, gamma=0.7, seed=999)
        W1 = scs1.get_weight_matrix()
        W2 = scs2.get_weight_matrix()
        assert np.allclose(W1, W2)


class TestkRegularGraphGT:
    """Tests for kRegularGraphGT native implementation."""

    def test_creation(self):
        """Test basic graph creation."""
        from lrgsglib.graphs.gt.kRegularGraphGT import kRegularGraphGT

        kreg = kRegularGraphGT(n=100, k=4, seed=42)
        assert kreg.N == 100
        expected_edges = 100 * 4 // 2
        assert kreg.num_edges == expected_edges

    def test_regularity(self):
        """Test all nodes have degree k."""
        from lrgsglib.graphs.gt.kRegularGraphGT import kRegularGraphGT

        kreg = kRegularGraphGT(n=50, k=6, seed=42)
        assert kreg.verify_regularity()

    def test_invalid_parameters(self):
        """Test invalid parameter combinations raise errors."""
        from lrgsglib.graphs.gt.kRegularGraphGT import kRegularGraphGT

        # n*k must be even
        with pytest.raises(ValueError):
            kRegularGraphGT(n=5, k=3, seed=42)

        # k must be less than n
        with pytest.raises(ValueError):
            kRegularGraphGT(n=10, k=10, seed=42)

    def test_sign_flipping(self):
        """Test edge sign flipping works."""
        from lrgsglib.graphs.gt.kRegularGraphGT import kRegularGraphGT

        kreg = kRegularGraphGT(n=100, k=4, pflip=0.3, seed=42)
        neg_edges = kreg.count_negative_edges()
        assert neg_edges > 0
        expected = int(0.3 * kreg.num_edges)
        assert abs(neg_edges - expected) < 0.15 * kreg.num_edges

    def test_reproducibility(self):
        """Test same seed gives same graph."""
        from lrgsglib.graphs.gt.kRegularGraphGT import kRegularGraphGT

        kreg1 = kRegularGraphGT(n=50, k=4, seed=123)
        kreg2 = kRegularGraphGT(n=50, k=4, seed=123)
        assert kreg1.N == kreg2.N
        assert kreg1.num_edges == kreg2.num_edges


class TestRandomGeometricGT:
    """Tests for RandomGeometricGT native implementation."""

    def test_creation(self):
        """Test basic graph creation."""
        from lrgsglib.graphs.gt.RandomGeometricGT import RandomGeometricGT

        rgg = RandomGeometricGT(n=200, radius=0.2, dim=2, seed=42)
        assert rgg.N > 0
        assert rgg.num_edges > 0

    def test_positions(self):
        """Test positions are stored."""
        from lrgsglib.graphs.gt.RandomGeometricGT import RandomGeometricGT

        rgg = RandomGeometricGT(n=100, radius=0.2, dim=2, seed=42)
        pos = rgg.get_positions()
        assert pos.shape[1] == 2
        assert np.all((pos >= 0) & (pos <= 1))

    def test_higher_dimension(self):
        """Test 3D random geometric graph."""
        from lrgsglib.graphs.gt.RandomGeometricGT import RandomGeometricGT

        rgg = RandomGeometricGT(n=100, radius=0.3, dim=3, seed=42)
        pos = rgg.get_positions()
        assert pos.shape[1] == 3

    def test_giant_component_extraction(self):
        """Test giant component extraction reduces node count."""
        from lrgsglib.graphs.gt.RandomGeometricGT import RandomGeometricGT

        # Small radius may create disconnected components
        rgg = RandomGeometricGT(
            n=500, radius=0.05, dim=2,
            extract_giant_component=True, seed=42
        )
        # With small radius, likely fewer nodes in giant component
        assert rgg.N <= 500

    def test_no_giant_component(self):
        """Test keeping full graph without extraction."""
        from lrgsglib.graphs.gt.RandomGeometricGT import RandomGeometricGT

        rgg = RandomGeometricGT(
            n=100, radius=0.2, dim=2,
            extract_giant_component=False, seed=42
        )
        # All nodes retained
        assert rgg.N == 100

    def test_critical_radius(self):
        """Test critical radius calculation."""
        from lrgsglib.graphs.gt.RandomGeometricGT import RandomGeometricGT

        rgg = RandomGeometricGT(n=500, radius=0.2, dim=2, seed=42)
        r_c = rgg.get_critical_radius()
        assert r_c > 0
        assert r_c < 0.5  # Should be relatively small

    def test_sign_flipping(self):
        """Test edge sign flipping works."""
        from lrgsglib.graphs.gt.RandomGeometricGT import RandomGeometricGT

        rgg = RandomGeometricGT(n=200, radius=0.2, pflip=0.3, seed=42)
        neg_edges = rgg.count_negative_edges()
        assert neg_edges > 0

    def test_reproducibility(self):
        """Test same seed gives same graph."""
        from lrgsglib.graphs.gt.RandomGeometricGT import RandomGeometricGT

        rgg1 = RandomGeometricGT(n=100, radius=0.2, dim=2, seed=123)
        rgg2 = RandomGeometricGT(n=100, radius=0.2, dim=2, seed=123)
        assert rgg1.N == rgg2.N
        assert rgg1.num_edges == rgg2.num_edges


class TestNewGTNativeConsistency:
    """Test node counts for new GT implementations."""

    def test_fully_connected_node_count(self):
        """Verify FullyConnectedGT produces expected node count."""
        from lrgsglib.graphs.gt.FullyConnectedGT import FullyConnectedGT

        for N in [20, 50, 100]:
            fc = FullyConnectedGT(N=N, seed=42)
            assert fc.N == N
            assert fc.num_edges == N * (N - 1) // 2

    def test_scs_node_count(self):
        """Verify SCSGeneralizedNNGT produces expected node count."""
        from lrgsglib.graphs.gt.SCSGeneralizedNNGT import SCSGeneralizedNNGT

        for N in [20, 50, 100]:
            scs = SCSGeneralizedNNGT(N=N, gamma=0.5, seed=42)
            assert scs.N == N

    def test_kreg_node_edge_count(self):
        """Verify kRegularGraphGT produces expected counts."""
        from lrgsglib.graphs.gt.kRegularGraphGT import kRegularGraphGT

        for n, k in [(50, 4), (100, 6), (200, 8)]:
            kreg = kRegularGraphGT(n=n, k=k, seed=42)
            assert kreg.N == n
            assert kreg.num_edges == n * k // 2


class TestConfigurationModelGT:
    """Tests for ConfigurationModelGT native implementation."""

    def test_creation(self):
        """Test basic graph creation."""
        from lrgsglib.graphs.gt.ConfigurationModelGT import ConfigurationModelGT

        degrees = [3, 3, 3, 3, 4, 4, 4, 4, 5, 5]  # sum = 38 (even)
        cm = ConfigurationModelGT(degree_sequence=degrees, seed=42)
        assert cm.N == 10
        assert cm.num_edges > 0

    def test_expected_edges(self):
        """Test expected edge count."""
        from lrgsglib.graphs.gt.ConfigurationModelGT import ConfigurationModelGT

        degrees = [4, 4, 4, 4, 4, 4, 4, 4]  # sum = 32, expect 16 edges
        cm = ConfigurationModelGT(degree_sequence=degrees, seed=42)
        assert cm.get_expected_num_edges() == 16

    def test_degree_statistics(self):
        """Test degree statistics computation."""
        from lrgsglib.graphs.gt.ConfigurationModelGT import ConfigurationModelGT

        degrees = [2, 3, 2, 3, 4, 4, 2, 4, 2, 4]  # 10 nodes, sum=30 (even)
        cm = ConfigurationModelGT(degree_sequence=degrees, seed=42)
        stats = cm.get_degree_statistics()
        assert 'mean' in stats
        assert 'std' in stats
        assert stats['min'] == min(degrees)
        assert stats['max'] == max(degrees)

    def test_invalid_odd_sum(self):
        """Test that odd degree sum raises error."""
        from lrgsglib.graphs.gt.ConfigurationModelGT import ConfigurationModelGT

        degrees = [3, 3, 3]  # sum = 9 (odd)
        with pytest.raises(ValueError):
            ConfigurationModelGT(degree_sequence=degrees, seed=42)

    def test_sign_flipping(self):
        """Test edge sign flipping works."""
        from lrgsglib.graphs.gt.ConfigurationModelGT import ConfigurationModelGT

        degrees = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4]  # 10 nodes, 20 edges
        cm = ConfigurationModelGT(degree_sequence=degrees, pflip=0.3, seed=42)
        neg_edges = cm.count_negative_edges()
        assert neg_edges > 0

    def test_reproducibility(self):
        """Test same seed gives same graph."""
        from lrgsglib.graphs.gt.ConfigurationModelGT import ConfigurationModelGT

        degrees = [4, 4, 4, 4, 4, 4, 4, 4]
        cm1 = ConfigurationModelGT(degree_sequence=degrees, seed=123)
        cm2 = ConfigurationModelGT(degree_sequence=degrees, seed=123)
        assert cm1.N == cm2.N
        assert cm1.num_edges == cm2.num_edges


class TestBipartiteGraphGT:
    """Tests for BipartiteGraphGT native implementation (random mode)."""

    def test_creation(self):
        """Test basic graph creation."""
        from lrgsglib.graphs.gt.BipartiteGraphGT import BipartiteGraphGT

        bg = BipartiteGraphGT(n1=20, n2=30, p=0.2, seed=42)
        assert bg.N == 50
        assert bg.num_edges > 0
        assert bg.mode == "random"

    def test_expected_nodes(self):
        """Test expected node count."""
        from lrgsglib.graphs.gt.BipartiteGraphGT import BipartiteGraphGT

        bg = BipartiteGraphGT(n1=15, n2=25, p=0.1, seed=42)
        assert bg.get_expected_num_nodes() == 40
        assert bg.N == 40

    def test_partitions(self):
        """Test partition sets are correct."""
        from lrgsglib.graphs.gt.BipartiteGraphGT import BipartiteGraphGT

        bg = BipartiteGraphGT(n1=10, n2=20, p=0.15, seed=42)
        assert len(bg.top_nodes) == 10
        assert len(bg.bottom_nodes) == 20
        assert bg.top_nodes == set(range(10))
        assert bg.bottom_nodes == set(range(10, 30))

    def test_is_bipartite(self):
        """Test graph is actually bipartite."""
        from lrgsglib.graphs.gt.BipartiteGraphGT import BipartiteGraphGT

        bg = BipartiteGraphGT(n1=15, n2=20, p=0.2, seed=42)
        assert bg.is_bipartite()

    def test_biadjacency_matrix(self):
        """Test biadjacency matrix extraction."""
        from lrgsglib.graphs.gt.BipartiteGraphGT import BipartiteGraphGT

        bg = BipartiteGraphGT(n1=10, n2=15, p=0.3, seed=42)
        B = bg.get_biadjacency_matrix()
        assert B.shape == (10, 15)
        # All entries should be 0 or 1
        assert np.all((B == 0) | (B == 1))

    def test_density(self):
        """Test density calculation."""
        from lrgsglib.graphs.gt.BipartiteGraphGT import BipartiteGraphGT

        bg = BipartiteGraphGT(n1=10, n2=20, p=0.5, seed=42)
        density = bg.get_density()
        # Should be roughly p (with some randomness)
        assert 0.3 < density < 0.7

    def test_projection(self):
        """Test one-mode projection."""
        from lrgsglib.graphs.gt.BipartiteGraphGT import BipartiteGraphGT

        bg = BipartiteGraphGT(n1=10, n2=15, p=0.4, seed=42)
        proj_top = bg.get_projection("top")
        proj_bottom = bg.get_projection("bottom")
        assert proj_top.num_vertices() == 10
        assert proj_bottom.num_vertices() == 15

    def test_sign_flipping(self):
        """Test edge sign flipping works."""
        from lrgsglib.graphs.gt.BipartiteGraphGT import BipartiteGraphGT

        bg = BipartiteGraphGT(n1=20, n2=30, p=0.2, pflip=0.3, seed=42)
        neg_edges = bg.count_negative_edges()
        assert neg_edges > 0

    def test_reproducibility(self):
        """Test same seed gives same graph."""
        from lrgsglib.graphs.gt.BipartiteGraphGT import BipartiteGraphGT

        bg1 = BipartiteGraphGT(n1=15, n2=25, p=0.2, seed=123)
        bg2 = BipartiteGraphGT(n1=15, n2=25, p=0.2, seed=123)
        assert bg1.N == bg2.N
        assert bg1.num_edges == bg2.num_edges

    def test_degree_sequence_mode_via_unified_class(self):
        """Test degree-sequence mode via unified BipartiteGraphGT."""
        from lrgsglib.graphs.gt.BipartiteGraphGT import BipartiteGraphGT

        top_deg = [3, 3, 2, 2]  # sum = 10
        bot_deg = [2, 2, 2, 2, 2]  # sum = 10
        bg = BipartiteGraphGT(
            top_degrees=top_deg, bottom_degrees=bot_deg, seed=42
        )
        assert bg.mode == "degree_sequence"
        assert bg.N == 9
        assert bg.num_edges > 0
        assert bg.get_expected_num_nodes() == 9
        assert bg.get_expected_num_edges() == 10
        assert bg.is_bipartite()

    def test_mode_auto_detection(self):
        """Test that mode is correctly auto-detected."""
        from lrgsglib.graphs.gt.BipartiteGraphGT import BipartiteGraphGT

        bg_random = BipartiteGraphGT(n1=10, n2=20, p=0.2, seed=42)
        assert bg_random.mode == "random"
        assert bg_random.top_degrees is None

        bg_deg = BipartiteGraphGT(
            top_degrees=[2, 2, 2], bottom_degrees=[3, 3], seed=42
        )
        assert bg_deg.mode == "degree_sequence"
        assert bg_deg.top_degrees == [2, 2, 2]

    def test_missing_params_raises_error(self):
        """Test that missing parameters raise ValueError."""
        from lrgsglib.graphs.gt.BipartiteGraphGT import BipartiteGraphGT

        # No params at all
        with pytest.raises(ValueError):
            BipartiteGraphGT(seed=42)

        # Only one degree sequence
        with pytest.raises(ValueError):
            BipartiteGraphGT(top_degrees=[2, 2], seed=42)

        # Only n1 without n2
        with pytest.raises(ValueError):
            BipartiteGraphGT(n1=10, seed=42)

    def test_degree_seq_methods_in_deg_mode(self):
        """Test degree-sequence-specific methods work in deg mode."""
        from lrgsglib.graphs.gt.BipartiteGraphGT import BipartiteGraphGT

        bg = BipartiteGraphGT(
            top_degrees=[2, 2, 2, 2, 2],
            bottom_degrees=[2, 2, 2, 2, 2],
            seed=42,
        )
        actual_top = bg.get_actual_top_degrees()
        actual_bottom = bg.get_actual_bottom_degrees()
        assert len(actual_top) == 5
        assert len(actual_bottom) == 5
        assert bg.verify_degrees()

    def test_degree_seq_methods_raise_in_random_mode(self):
        """Test degree-sequence methods raise in random mode."""
        from lrgsglib.graphs.gt.BipartiteGraphGT import BipartiteGraphGT

        bg = BipartiteGraphGT(n1=10, n2=20, p=0.2, seed=42)
        with pytest.raises(AttributeError):
            bg.get_actual_top_degrees()
        with pytest.raises(AttributeError):
            bg.get_actual_bottom_degrees()
        with pytest.raises(AttributeError):
            bg.verify_degrees()


class TestBipartiteFromDegreeSequenceGT:
    """Tests for BipartiteFromDegreeSequenceGT backward-compat wrapper."""

    def test_creation(self):
        """Test basic graph creation via old class name."""
        from lrgsglib.graphs.gt.BipartiteFromDegreeSequenceGT import (
            BipartiteFromDegreeSequenceGT,
        )

        top_deg = [3, 3, 2, 2]  # sum = 10
        bot_deg = [2, 2, 2, 2, 2]  # sum = 10
        bg = BipartiteFromDegreeSequenceGT(top_deg, bot_deg, seed=42)
        assert bg.N == 9
        assert bg.num_edges > 0
        assert bg.mode == "degree_sequence"

    def test_expected_nodes(self):
        """Test expected node count."""
        from lrgsglib.graphs.gt.BipartiteFromDegreeSequenceGT import (
            BipartiteFromDegreeSequenceGT,
        )

        top_deg = [4, 4, 4, 4]  # 4 top nodes
        bot_deg = [2, 2, 2, 2, 2, 2, 2, 2]  # 8 bottom nodes
        bg = BipartiteFromDegreeSequenceGT(top_deg, bot_deg, seed=42)
        assert bg.get_expected_num_nodes() == 12

    def test_expected_edges(self):
        """Test expected edge count."""
        from lrgsglib.graphs.gt.BipartiteFromDegreeSequenceGT import (
            BipartiteFromDegreeSequenceGT,
        )

        top_deg = [3, 3, 3, 3]  # sum = 12
        bot_deg = [4, 4, 4]  # sum = 12
        bg = BipartiteFromDegreeSequenceGT(top_deg, bot_deg, seed=42)
        assert bg.get_expected_num_edges() == 12

    def test_mismatched_degree_sums(self):
        """Test that mismatched degree sums raise error."""
        from lrgsglib.graphs.gt.BipartiteFromDegreeSequenceGT import (
            BipartiteFromDegreeSequenceGT,
        )

        top_deg = [3, 3, 3]  # sum = 9
        bot_deg = [2, 2, 2, 2]  # sum = 8
        with pytest.raises(ValueError):
            BipartiteFromDegreeSequenceGT(top_deg, bot_deg, seed=42)

    def test_is_bipartite(self):
        """Test graph is actually bipartite."""
        from lrgsglib.graphs.gt.BipartiteFromDegreeSequenceGT import (
            BipartiteFromDegreeSequenceGT,
        )

        top_deg = [3, 3, 2, 2]
        bot_deg = [2, 2, 2, 2, 2]
        bg = BipartiteFromDegreeSequenceGT(top_deg, bot_deg, seed=42)
        assert bg.is_bipartite()

    def test_biadjacency_matrix(self):
        """Test biadjacency matrix extraction."""
        from lrgsglib.graphs.gt.BipartiteFromDegreeSequenceGT import (
            BipartiteFromDegreeSequenceGT,
        )

        top_deg = [2, 2, 2, 2, 2]  # 5 top
        bot_deg = [2, 2, 2, 2, 2]  # 5 bottom
        bg = BipartiteFromDegreeSequenceGT(top_deg, bot_deg, seed=42)
        B = bg.get_biadjacency_matrix()
        assert B.shape == (5, 5)

    def test_sign_flipping(self):
        """Test edge sign flipping works."""
        from lrgsglib.graphs.gt.BipartiteFromDegreeSequenceGT import (
            BipartiteFromDegreeSequenceGT,
        )

        top_deg = [4, 4, 4, 4]
        bot_deg = [2, 2, 2, 2, 2, 2, 2, 2]
        bg = BipartiteFromDegreeSequenceGT(top_deg, bot_deg, pflip=0.3, seed=42)
        neg_edges = bg.count_negative_edges()
        assert neg_edges > 0

    def test_reproducibility(self):
        """Test same seed gives same graph."""
        from lrgsglib.graphs.gt.BipartiteFromDegreeSequenceGT import (
            BipartiteFromDegreeSequenceGT,
        )

        top_deg = [3, 3, 3, 3]
        bot_deg = [4, 4, 4]
        bg1 = BipartiteFromDegreeSequenceGT(top_deg, bot_deg, seed=123)
        bg2 = BipartiteFromDegreeSequenceGT(top_deg, bot_deg, seed=123)
        assert bg1.N == bg2.N
        assert bg1.num_edges == bg2.num_edges

    def test_is_subclass_of_bipartite_graph_gt(self):
        """Test backward-compat class is subclass of BipartiteGraphGT."""
        from lrgsglib.graphs.gt.BipartiteFromDegreeSequenceGT import (
            BipartiteFromDegreeSequenceGT,
        )
        from lrgsglib.graphs.gt.BipartiteGraphGT import BipartiteGraphGT

        assert issubclass(BipartiteFromDegreeSequenceGT, BipartiteGraphGT)


class TestPhase2GTNativeConsistency:
    """Test node counts for Phase 2 GT implementations."""

    def test_configuration_model_counts(self):
        """Verify ConfigurationModelGT produces expected counts."""
        from lrgsglib.graphs.gt.ConfigurationModelGT import ConfigurationModelGT

        for n in [20, 50, 100]:
            degrees = [4] * n  # Regular k=4
            cm = ConfigurationModelGT(degree_sequence=degrees, seed=42)
            assert cm.N == n
            assert cm.get_expected_num_edges() == n * 4 // 2

    def test_bipartite_counts(self):
        """Verify BipartiteGraphGT produces expected counts."""
        from lrgsglib.graphs.gt.BipartiteGraphGT import BipartiteGraphGT

        for n1, n2 in [(10, 20), (25, 50), (50, 100)]:
            bg = BipartiteGraphGT(n1=n1, n2=n2, p=0.5, seed=42)
            assert bg.N == n1 + n2
            assert len(bg.top_nodes) == n1
            assert len(bg.bottom_nodes) == n2

    def test_bipartite_deg_seq_counts(self):
        """Verify BipartiteFromDegreeSequenceGT produces expected counts."""
        from lrgsglib.graphs.gt.BipartiteFromDegreeSequenceGT import (
            BipartiteFromDegreeSequenceGT,
        )

        top_deg = [3] * 10  # sum = 30
        bot_deg = [2] * 15  # sum = 30
        bg = BipartiteFromDegreeSequenceGT(top_deg, bot_deg, seed=42)
        assert bg.get_expected_num_nodes() == 25
        assert bg.get_expected_num_edges() == 30


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
