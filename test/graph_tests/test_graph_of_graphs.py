"""
Tests for GraphOfGraphs class.

Tests the generalized hierarchical graph structure implementation for both
NetworkX and graph-tool backends.
"""

import numpy as np
import pytest

# Skip all tests if graph-tool not available
gt_available = True
try:
    import graph_tool as gt
except ImportError:
    gt_available = False


class TestGraphOfGraphsNX:
    """Test NetworkX implementation of GraphOfGraphs."""

    def test_basic_creation(self):
        """Test basic graph creation with Lattice2D base and ErdosRenyi fibers."""
        from lrgsglib.graphs.nx.GraphOfGraphsNX import GraphOfGraphsNX

        gog = GraphOfGraphsNX(
            base_graph_type="Lattice2D",
            base_params={"side1": 4, "geo": "sqr"},
            fiber_graph_type="ErdosRenyi",
            fiber_params={"n": 5, "p": 0.5},
            seed=42,
        )

        # Check dimensions
        assert gog.N_base == 16  # 4x4 lattice
        assert gog.N_fiber == 5
        assert gog.N == 16 + 16 * 5  # 96 total nodes

    def test_metadata_structure(self):
        """Test that gog_structure metadata is properly populated."""
        from lrgsglib.graphs.nx.GraphOfGraphsNX import GraphOfGraphsNX

        # Use Lattice2D for fiber to get predictable size
        gog = GraphOfGraphsNX(
            base_graph_type="Lattice2D",
            base_params={"side1": 3},
            fiber_graph_type="Lattice2D",
            fiber_params={"side1": 2},  # 2x2 = 4 nodes
            anchor_policy="first",
            seed=123,
        )

        meta = gog.gog_structure
        assert meta["base_graph_class"] == "Lattice2D"
        assert meta["fiber_graph_class"] == "Lattice2D"
        assert meta["N_base"] == 9
        assert meta["N_fiber"] == 4
        assert meta["N_total"] == 9 + 9 * 4
        assert meta["anchor_policy"] == "first"
        assert not meta["heterogeneous"]
        assert meta["can_use_separated_spectrum"] is True

    def test_vertex_indices(self):
        """Test vertex index methods."""
        from lrgsglib.graphs.nx.GraphOfGraphsNX import GraphOfGraphsNX

        gog = GraphOfGraphsNX(
            base_graph_type="Lattice2D",
            base_params={"side1": 2},
            fiber_graph_type="ErdosRenyi",
            fiber_params={"n": 3, "p": 0.5},
            seed=42,
        )

        # Base vertices should be [0, 1, 2, 3]
        assert gog.base_vertex_indices() == [0, 1, 2, 3]

        # Fiber vertices for base 0: [4, 5, 6]
        assert gog.fiber_vertex_indices(0) == [4, 5, 6]

        # Fiber vertices for base 1: [7, 8, 9]
        assert gog.fiber_vertex_indices(1) == [7, 8, 9]

        # Anchor for base 0 with 'first' policy
        assert gog.anchor_vertex(0) == 4

    def test_vertex_layer(self):
        """Test vertex_layer method."""
        from lrgsglib.graphs.nx.GraphOfGraphsNX import GraphOfGraphsNX

        gog = GraphOfGraphsNX(
            base_graph_type="Lattice2D",
            base_params={"side1": 2},
            fiber_graph_type="ErdosRenyi",
            fiber_params={"n": 3, "p": 0.5},
            seed=42,
        )

        # First 4 are base
        for v in range(4):
            assert gog.vertex_layer(v) == "base"

        # Rest are fiber
        for v in range(4, gog.N):
            assert gog.vertex_layer(v) == "fiber"

    def test_anchor_policies(self):
        """Test different anchor policies."""
        from lrgsglib.graphs.nx.GraphOfGraphsNX import GraphOfGraphsNX

        n_fiber = 5
        base_params = {"side1": 2}
        fiber_params = {"n": n_fiber, "p": 0.5}

        # First policy
        gog_first = GraphOfGraphsNX(
            base_graph_type="Lattice2D",
            base_params=base_params,
            fiber_graph_type="ErdosRenyi",
            fiber_params=fiber_params,
            anchor_policy="first",
            seed=42,
        )
        # Anchor should be first node of fiber (index 0)
        n_base = gog_first.N_base
        for i in range(n_base):
            assert gog_first.anchor_vertex(i) == n_base + i * n_fiber + 0

        # Center policy
        gog_center = GraphOfGraphsNX(
            base_graph_type="Lattice2D",
            base_params=base_params,
            fiber_graph_type="ErdosRenyi",
            fiber_params=fiber_params,
            anchor_policy="center",
            seed=42,
        )
        for i in range(n_base):
            assert gog_center.anchor_vertex(i) == n_base + i * n_fiber + n_fiber // 2

        # Last policy
        gog_last = GraphOfGraphsNX(
            base_graph_type="Lattice2D",
            base_params=base_params,
            fiber_graph_type="ErdosRenyi",
            fiber_params=fiber_params,
            anchor_policy="last",
            seed=42,
        )
        for i in range(n_base):
            assert gog_last.anchor_vertex(i) == n_base + i * n_fiber + n_fiber - 1

    def test_callable_anchor_policy(self):
        """Test custom callable anchor policy."""
        from lrgsglib.graphs.nx.GraphOfGraphsNX import GraphOfGraphsNX

        # Custom policy: alternate between first and last
        def alternating_anchor(base_idx, n_fiber):
            return 0 if base_idx % 2 == 0 else n_fiber - 1

        gog = GraphOfGraphsNX(
            base_graph_type="Lattice2D",
            base_params={"side1": 2},
            fiber_graph_type="ErdosRenyi",
            fiber_params={"n": 5, "p": 0.5},
            anchor_policy=alternating_anchor,
            seed=42,
        )

        n_base = gog.N_base
        n_fiber = gog.N_fiber

        # Check alternating pattern
        for i in range(n_base):
            expected_local = 0 if i % 2 == 0 else n_fiber - 1
            expected_global = n_base + i * n_fiber + expected_local
            assert gog.anchor_vertex(i) == expected_global

    def test_can_use_separated_spectrum_homogeneous(self):
        """Test separated spectrum flag for homogeneous fibers."""
        from lrgsglib.graphs.nx.GraphOfGraphsNX import GraphOfGraphsNX

        gog = GraphOfGraphsNX(
            base_graph_type="Lattice2D",
            base_params={"side1": 3},
            fiber_graph_type="ErdosRenyi",
            fiber_params={"n": 5, "p": 0.3},
            anchor_policy="first",
            seed=42,
        )

        assert gog.can_use_separated_spectrum() is True

    def test_can_use_separated_spectrum_random_anchor(self):
        """Test that random anchor policy disables separated spectrum."""
        from lrgsglib.graphs.nx.GraphOfGraphsNX import GraphOfGraphsNX

        gog = GraphOfGraphsNX(
            base_graph_type="Lattice2D",
            base_params={"side1": 3},
            fiber_graph_type="ErdosRenyi",
            fiber_params={"n": 5, "p": 0.3},
            anchor_policy="random",
            seed=42,
        )

        # Random anchors are not uniform, so separated spectrum is invalid
        # (Unless by chance all random values are the same)
        anchor_set = set(gog._anchor_indices)
        if len(anchor_set) > 1:
            assert gog.can_use_separated_spectrum() is False

    def test_dirac_structure_compatibility(self):
        """Test dirac_structure property for backward compatibility."""
        from lrgsglib.graphs.nx.GraphOfGraphsNX import GraphOfGraphsNX

        # Use Lattice2D for fiber to get predictable size
        gog = GraphOfGraphsNX(
            base_graph_type="Lattice2D",
            base_params={"side1": 3},
            fiber_graph_type="Lattice2D",
            fiber_params={"side1": 2},
            seed=42,
        )

        dirac = gog.dirac_structure
        assert dirac["base_nodes"] == 9
        assert dirac["fiber_nodes"] == 4
        assert dirac["total_nodes"] == 9 + 9 * 4

    def test_edge_count(self):
        """Test that edge count is correct."""
        from lrgsglib.graphs.nx.GraphOfGraphsNX import GraphOfGraphsNX

        # Use fully connected fibers for predictable edge count
        gog = GraphOfGraphsNX(
            base_graph_type="Lattice2D",
            base_params={"side1": 2, "geo": "sqr", "pbc": False},
            fiber_graph_type="ErdosRenyi",
            fiber_params={"n": 3, "p": 1.0},  # Complete graph
            seed=42,
        )

        n_base = gog.N_base
        n_fiber = gog.N_fiber

        # Base edges: 2x2 non-periodic square = 4 edges
        base_edges = 4

        # Fiber edges per fiber: complete graph on 3 nodes = 3 edges
        # Total fiber edges: 4 fibers * 3 edges = 12
        fiber_edges = n_base * 3

        # Anchor edges: 4 (one per base node)
        anchor_edges = n_base

        expected_total = base_edges + fiber_edges + anchor_edges
        assert gog.G.number_of_edges() == expected_total


class TestGraphOfGraphsGT:
    """Test graph-tool implementation of GraphOfGraphs."""

    @pytest.mark.skipif(not gt_available, reason="graph-tool not installed")
    def test_basic_creation(self):
        """Test basic graph creation with GT backend."""
        from lrgsglib.graphs.gt.GraphOfGraphsGT import GraphOfGraphsGT

        # Use Lattice2D for fiber to get predictable size
        gog = GraphOfGraphsGT(
            base_graph_type="Lattice2D",
            base_params={"side1": 4, "geo": "sqr"},
            fiber_graph_type="Lattice2D",
            fiber_params={"side1": 2},  # 2x2 = 4 nodes
            seed=42,
        )

        # Check dimensions
        assert gog.N_base == 16  # 4x4 lattice
        assert gog.N_fiber == 4
        assert gog.N == 16 + 16 * 4  # 80 total nodes

    @pytest.mark.skipif(not gt_available, reason="graph-tool not installed")
    def test_metadata_structure(self):
        """Test that gog_structure metadata is properly populated."""
        from lrgsglib.graphs.gt.GraphOfGraphsGT import GraphOfGraphsGT

        # Use Lattice2D for fiber to get predictable size
        gog = GraphOfGraphsGT(
            base_graph_type="Lattice2D",
            base_params={"side1": 3},
            fiber_graph_type="Lattice2D",
            fiber_params={"side1": 2},
            anchor_policy="first",
            seed=123,
        )

        meta = gog.gog_structure
        assert meta["base_graph_class"] == "Lattice2D"
        assert meta["fiber_graph_class"] == "Lattice2D"
        assert meta["N_base"] == 9
        assert meta["N_fiber"] == 4
        assert meta["N_total"] == 9 + 9 * 4

    @pytest.mark.skipif(not gt_available, reason="graph-tool not installed")
    def test_vertex_layer(self):
        """Test vertex_layer method for GT backend."""
        from lrgsglib.graphs.gt.GraphOfGraphsGT import GraphOfGraphsGT

        gog = GraphOfGraphsGT(
            base_graph_type="Lattice2D",
            base_params={"side1": 2},
            fiber_graph_type="ErdosRenyi",
            fiber_params={"n": 3, "p": 0.5},
            seed=42,
        )

        # First 4 are base
        for v in range(4):
            assert gog.vertex_layer(v) == "base"

        # Rest are fiber
        for v in range(4, gog.N):
            assert gog.vertex_layer(v) == "fiber"


class TestGraphOfGraphsFacade:
    """Test the unified GraphOfGraphs facade."""

    def test_nx_engine(self):
        """Test facade with NX engine."""
        from lrgsglib.graphs import GraphOfGraphs

        # Use Lattice2D for fiber to get predictable size
        gog = GraphOfGraphs(
            base_graph_type="Lattice2D",
            base_params={"side1": 3},
            fiber_graph_type="Lattice2D",
            fiber_params={"side1": 2},
            engine="nx",
            seed=42,
        )

        assert gog.N_base == 9
        assert gog.N_fiber == 4

    @pytest.mark.skipif(not gt_available, reason="graph-tool not installed")
    def test_gt_engine(self):
        """Test facade with GT engine."""
        from lrgsglib.graphs import GraphOfGraphs

        # Use Lattice2D for fiber to get predictable size
        gog = GraphOfGraphs(
            base_graph_type="Lattice2D",
            base_params={"side1": 3},
            fiber_graph_type="Lattice2D",
            fiber_params={"side1": 2},
            engine="gt",
            seed=42,
        )

        assert gog.N_base == 9
        assert gog.N_fiber == 4

    @pytest.mark.skipif(not gt_available, reason="graph-tool not installed")
    def test_nx_gt_equivalence(self):
        """Test that NX and GT produce equivalent structures."""
        from lrgsglib.graphs import GraphOfGraphs

        # Use Lattice2D for predictable sizes
        gog_nx = GraphOfGraphs(
            base_graph_type="Lattice2D",
            base_params={"side1": 3},
            fiber_graph_type="Lattice2D",
            fiber_params={"side1": 2},
            engine="nx",
            seed=42,
        )

        gog_gt = GraphOfGraphs(
            base_graph_type="Lattice2D",
            base_params={"side1": 3},
            fiber_graph_type="Lattice2D",
            fiber_params={"side1": 2},
            engine="gt",
            seed=42,
        )

        # Dimensions should match exactly
        assert gog_nx.N == gog_gt.N
        assert gog_nx.N_base == gog_gt.N_base
        assert gog_nx.N_fiber == gog_gt.N_fiber


class TestSpectralMethods:
    """Test spectral computation methods for GraphOfGraphs."""

    def test_separated_spectrum_basic(self):
        """Test basic separated spectrum computation."""
        from lrgsglib.graphs.nx.GraphOfGraphsNX import GraphOfGraphsNX

        gog = GraphOfGraphsNX(
            base_graph_type="Lattice2D",
            base_params={"side1": 3, "pbc": False},
            fiber_graph_type="Lattice2D",
            fiber_params={"side1": 2, "pbc": False},
            anchor_policy="first",
            seed=42,
        )

        assert gog.can_use_separated_spectrum()

        # Compute separated spectrum
        spectrum = gog.compute_separated_spectrum()

        # The separated spectrum formula gives N_base * N_fiber eigenvalues
        # (not N_total = N_base + N_base * N_fiber)
        # This is an efficient approximation that captures spectral structure
        assert len(spectrum) == gog.N_base * gog.N_fiber

        # First eigenvalue should be ~0
        assert abs(spectrum[0]) < 1e-10

        # All eigenvalues should be non-negative (Laplacian property)
        assert all(ev >= -1e-10 for ev in spectrum)

    def test_separated_spectrum_count(self):
        """Test that separated spectrum has correct count (N_base * N_fiber)."""
        from lrgsglib.graphs.nx.GraphOfGraphsNX import GraphOfGraphsNX

        gog = GraphOfGraphsNX(
            base_graph_type="Lattice2D",
            base_params={"side1": 2, "pbc": False},
            fiber_graph_type="Lattice2D",
            fiber_params={"side1": 3, "pbc": False},
            anchor_policy="first",
            seed=42,
        )

        spectrum = gog.compute_separated_spectrum()

        # 2x2 base = 4 nodes, 3x3 fiber = 9 nodes
        # Separated spectrum: 4 * 9 = 36 eigenvalues
        assert len(spectrum) == 4 * 9

    def test_base_fiber_dimensions(self):
        """Test get_base_fiber_dimensions method."""
        from lrgsglib.graphs.nx.GraphOfGraphsNX import GraphOfGraphsNX

        # Use Lattice2D for fiber to get predictable size
        gog = GraphOfGraphsNX(
            base_graph_type="Lattice2D",
            base_params={"side1": 3},
            fiber_graph_type="Lattice2D",
            fiber_params={"side1": 2},  # 2x2 = 4 nodes
            seed=42,
        )

        n_base, n_fiber = gog.get_base_fiber_dimensions()
        assert n_base == 9
        assert n_fiber == 4


class TestDiracEquivalence:
    """Test equivalence with existing DiracComb/DiracBrush structures."""

    def test_dirac_comb_equivalence(self):
        """Test that GraphOfGraphs with 1D base + 1D fiber matches DiracComb."""
        from lrgsglib.graphs.nx.GraphOfGraphsNX import GraphOfGraphsNX
        from lrgsglib.graphs.nx.DiracLatticeNX import DiracCombGraphNX

        base_nodes = 5
        fiber_nodes = 4

        # Create DiracComb
        dc = DiracCombGraphNX(
            base_nodes=base_nodes,
            fiber_nodes=fiber_nodes,
            periodic=True,
            seed=42,
        )

        # Create equivalent GraphOfGraphs
        # Note: DiracComb uses line/ring, which is Lattice2D with side2=1
        gog = GraphOfGraphsNX(
            base_graph_type="Lattice2D",
            base_params={"side1": base_nodes, "side2": 1, "pbc": True},
            fiber_graph_type="Lattice2D",
            fiber_params={"side1": fiber_nodes, "side2": 1, "pbc": True},
            anchor_policy="first",
            seed=42,
        )

        # Dimensions should match
        assert gog.N == dc.N
        assert gog.N_base == dc.base_nodes
        assert gog.N_fiber == dc.fiber_nodes

        # Edge counts should match
        assert gog.G.number_of_edges() == dc.G.number_of_edges()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
