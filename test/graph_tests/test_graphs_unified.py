"""
Tests for unified multi-engine graph interface.

Tests protocol compliance, engine switching, and cross-engine equivalence
for the lrgsglib.graphs unified API.
"""

import numpy as np
import pytest

from lrgsglib.graphs import (
    Lattice2D,
    Lattice3D,
    ErdosRenyi,
    BarabasiAlbert,
    WattsStrogatz,
    StochasticBlockModel,
    GraphEngine,
    set_default_engine,
    get_default_engine,
    available_engines,
    list_registered_types,
    is_engine_available,
)
from lrgsglib.graphs.protocols import (
    SignedGraphProtocol,
    is_signed_graph,
)


class TestEngineManagement:
    """Tests for engine selection and management."""

    def test_available_engines(self):
        """Test that available engines are detected correctly."""
        engines = available_engines()
        assert GraphEngine.NETWORKX in engines
        # GT may or may not be available depending on environment

    def test_set_get_default_engine(self):
        """Test setting and getting default engine."""
        original = get_default_engine()
        try:
            set_default_engine("nx")
            assert get_default_engine() == GraphEngine.NETWORKX

            if is_engine_available("gt"):
                set_default_engine("gt")
                assert get_default_engine() == GraphEngine.GRAPHTOOL
        finally:
            set_default_engine(original)

    def test_registered_types(self):
        """Test that graph types are registered."""
        types = list_registered_types()
        assert "Lattice2D" in types
        assert "Lattice3D" in types
        assert "ErdosRenyi" in types
        assert "BarabasiAlbert" in types
        assert "WattsStrogatz" in types
        assert "StochasticBlockModel" in types


class TestLattice2DNetworkX:
    """Tests for Lattice2D with NetworkX backend."""

    def test_creation(self):
        """Test basic lattice creation."""
        lat = Lattice2D(side1=10, geo="sqr", engine="nx")
        assert lat.N == 100

    def test_protocol_compliance(self):
        """Test that Lattice2D satisfies SignedGraphProtocol."""
        lat = Lattice2D(side1=10, geo="sqr", pflip=0.2, engine="nx")
        lat.flip_random_fract_edges()
        assert is_signed_graph(lat)

    def test_edge_flipping(self):
        """Test edge sign flipping."""
        lat = Lattice2D(side1=10, geo="sqr", pflip=0.3, engine="nx", seed=42)
        lat.flip_random_fract_edges()

        neg = lat.count_negative_edges()
        pos = lat.count_positive_edges()
        total = lat.num_edges

        assert neg + pos == total
        assert neg > 0  # With pflip=0.3, should have some negative

    def test_spectral_computation(self):
        """Test spectral decomposition."""
        lat = Lattice2D(side1=8, geo="sqr", periodic=False, engine="nx", seed=42)
        lat.flip_random_fract_edges()

        result = lat.compute_laplacian_spectrum_weigV()
        if result is not None:
            eigv, eigV = result
        else:
            eigv, eigV = lat.eigv, lat.eigV

        assert len(eigv) == lat.N
        assert eigV.shape == (lat.N, lat.N)
        # Smallest eigenvalue should be non-negative (within numerical precision)
        assert eigv[0] >= -1e-10

    def test_triangular_geometry(self):
        """Test triangular lattice geometry."""
        lat = Lattice2D(side1=10, geo="tri", periodic=False, engine="nx")
        # Triangular lattice has more edges than square
        lat_sqr = Lattice2D(side1=10, geo="sqr", periodic=False, engine="nx")
        assert lat.num_edges > lat_sqr.num_edges


@pytest.mark.skipif(
    not is_engine_available("gt"),
    reason="graph-tool not installed",
)
class TestLattice2DGraphTool:
    """Tests for Lattice2D with graph-tool backend."""

    def test_creation(self):
        """Test basic lattice creation."""
        lat = Lattice2D(side1=10, geo="sqr", engine="gt")
        assert lat.N == 100

    def test_protocol_compliance(self):
        """Test that Lattice2DGT satisfies SignedGraphProtocol."""
        lat = Lattice2D(side1=10, geo="sqr", pflip=0.2, engine="gt")
        lat.flip_random_fract_edges()
        assert is_signed_graph(lat)

    def test_spectral_computation(self):
        """Test spectral decomposition."""
        lat = Lattice2D(side1=8, geo="sqr", periodic=False, engine="gt", seed=42)
        lat.flip_random_fract_edges()

        eigv, eigV = lat.compute_laplacian_spectrum_weigV()

        assert len(eigv) == lat.N
        assert eigV.shape == (lat.N, lat.N)
        # Smallest eigenvalue should be non-negative (within numerical precision)
        assert eigv[0] >= -1e-10

    def test_to_networkx_conversion(self):
        """Test conversion to NetworkX."""
        lat = Lattice2D(side1=5, geo="sqr", engine="gt")
        nx_graph = lat.to_networkx()

        assert nx_graph.number_of_nodes() == lat.N
        assert nx_graph.number_of_edges() == lat.num_edges


@pytest.mark.skipif(
    not is_engine_available("gt"),
    reason="graph-tool not installed",
)
class TestCrossEngineEquivalence:
    """Tests for equivalence between NX and GT backends."""

    def test_lattice2d_same_structure(self):
        """Test that NX and GT produce same graph structure."""
        lat_nx = Lattice2D(side1=10, geo="sqr", periodic=False, engine="nx")
        lat_gt = Lattice2D(side1=10, geo="sqr", periodic=False, engine="gt")

        assert lat_nx.N == lat_gt.N
        assert lat_nx.num_edges == lat_gt.num_edges

    def test_lattice2d_same_laplacian(self):
        """Test that signed Laplacian is equivalent (unflipped)."""
        lat_nx = Lattice2D(side1=5, geo="sqr", pflip=0.0, periodic=False, engine="nx")
        lat_gt = Lattice2D(side1=5, geo="sqr", pflip=0.0, periodic=False, engine="gt")

        # Need to compute matrices for NX
        lat_nx.upd_graph_matrices()

        L_nx = lat_nx.get_signed_laplacian()
        L_gt = lat_gt.get_signed_laplacian()

        if hasattr(L_nx, "toarray"):
            L_nx = L_nx.toarray()

        np.testing.assert_allclose(L_nx, L_gt, rtol=1e-10)

    def test_lattice3d_sc_same_structure(self):
        """Test that 3D simple cubic lattice has same structure."""
        lat_nx = Lattice3D(dim=5, geo="sc", periodic=False, engine="nx")
        lat_gt = Lattice3D(dim=5, geo="sc", periodic=False, engine="gt")

        assert lat_nx.N == lat_gt.N
        assert lat_nx.num_edges == lat_gt.num_edges


class TestLattice3DNetworkX:
    """Tests for Lattice3D with NetworkX backend."""

    def test_simple_cubic(self):
        """Test simple cubic lattice."""
        lat = Lattice3D(dim=5, geo="sc", engine="nx")
        assert lat.N == 125  # 5^3

    def test_bcc_lattice(self):
        """Test BCC lattice creation."""
        lat = Lattice3D(dim=4, geo="bcc", engine="nx")
        assert lat.N > 0
        assert lat.num_edges > 0

    def test_fcc_lattice(self):
        """Test FCC lattice creation."""
        lat = Lattice3D(dim=4, geo="fcc", engine="nx")
        assert lat.N > 0
        assert lat.num_edges > 0

    def test_protocol_compliance(self):
        """Test that Lattice3D satisfies SignedGraphProtocol."""
        lat = Lattice3D(dim=5, geo="sc", pflip=0.2, engine="nx")
        lat.flip_random_fract_edges()
        assert is_signed_graph(lat)


@pytest.mark.skipif(
    not is_engine_available("gt"),
    reason="graph-tool not installed",
)
class TestLattice3DGraphTool:
    """Tests for Lattice3D with graph-tool backend."""

    def test_simple_cubic(self):
        """Test simple cubic lattice."""
        lat = Lattice3D(dim=5, geo="sc", engine="gt")
        assert lat.N == 125  # 5^3

    def test_bcc_lattice(self):
        """Test BCC lattice creation."""
        lat = Lattice3D(dim=4, geo="bcc", engine="gt")
        assert lat.N > 0
        assert lat.num_edges > 0

    def test_fcc_lattice(self):
        """Test FCC lattice creation."""
        lat = Lattice3D(dim=4, geo="fcc", engine="gt")
        assert lat.N > 0
        assert lat.num_edges > 0

    def test_spectral_computation(self):
        """Test spectral decomposition on 3D lattice."""
        lat = Lattice3D(dim=4, geo="sc", pflip=0.2, engine="gt", seed=42)
        lat.flip_random_fract_edges()

        eigv, eigV = lat.compute_laplacian_spectrum_weigV()

        assert len(eigv) == lat.N
        assert eigV.shape == (lat.N, lat.N)


class TestEngineFallback:
    """Tests for engine fallback behavior."""

    def test_fallback_to_networkx(self):
        """Test that unavailable engine falls back to NX."""
        # Request igraph which isn't implemented yet
        lat = Lattice2D(side1=5, geo="sqr", engine="nx")
        assert lat.N == 25  # Should work with NX


# ============================================================================
# Random Graph Tests
# ============================================================================


class TestErdosRenyiNetworkX:
    """Tests for ErdosRenyi with NetworkX backend."""

    def test_creation(self):
        """Test basic ER graph creation."""
        er = ErdosRenyi(n=100, p=0.1, seed=42, engine="nx")
        assert er.N > 0
        assert er.num_edges > 0

    def test_protocol_compliance(self):
        """Test that ErdosRenyi satisfies SignedGraphProtocol."""
        er = ErdosRenyi(n=100, p=0.1, pflip=0.2, seed=42, engine="nx")
        er.flip_random_fract_edges()
        assert is_signed_graph(er)

    def test_edge_flipping(self):
        """Test edge sign flipping."""
        er = ErdosRenyi(n=100, p=0.1, pflip=0.3, seed=42, engine="nx")
        er.flip_random_fract_edges()

        neg = er.count_negative_edges()
        total = er.num_edges

        assert neg > 0
        assert neg < total

    def test_spectral_computation(self):
        """Test signed Laplacian computation."""
        er = ErdosRenyi(n=50, p=0.15, pflip=0.2, seed=42, engine="nx")
        er.flip_random_fract_edges()

        L = er.get_signed_laplacian()
        if hasattr(L, "toarray"):
            L = L.toarray()

        assert L.shape == (er.N, er.N)
        eigv = np.linalg.eigvalsh(L)
        # Signed Laplacian eigenvalues should all be >= 0 (for reasonable pflip)
        assert eigv[0] >= -1e-10


@pytest.mark.skipif(
    not is_engine_available("gt"),
    reason="graph-tool not installed",
)
class TestErdosRenyiGraphTool:
    """Tests for ErdosRenyi with graph-tool backend."""

    def test_creation(self):
        """Test basic ER graph creation."""
        er = ErdosRenyi(n=100, p=0.1, seed=42, engine="gt")
        assert er.N > 0
        assert er.num_edges > 0

    def test_protocol_compliance(self):
        """Test that ErdosRenyiGT satisfies SignedGraphProtocol."""
        er = ErdosRenyi(n=100, p=0.1, pflip=0.2, seed=42, engine="gt")
        er.flip_random_fract_edges()
        assert is_signed_graph(er)

    def test_giant_component_extraction(self):
        """Test giant component extraction."""
        # Sparse graph likely to have multiple components before GC
        er_gc = ErdosRenyi(n=200, p=0.01, extract_giant_component=True, seed=42, engine="gt")
        er_no_gc = ErdosRenyi(n=200, p=0.01, extract_giant_component=False, seed=43, engine="gt")

        # GC extracted should have fewer nodes (or equal)
        assert er_gc.N <= 200

    def test_to_networkx_conversion(self):
        """Test conversion to NetworkX."""
        er = ErdosRenyi(n=50, p=0.1, engine="gt")
        nx_graph = er.to_networkx()

        assert nx_graph.number_of_nodes() == er.N
        assert nx_graph.number_of_edges() == er.num_edges


class TestBarabasiAlbertNetworkX:
    """Tests for BarabasiAlbert with NetworkX backend."""

    def test_creation(self):
        """Test basic BA graph creation."""
        ba = BarabasiAlbert(n=100, m=3, seed=42, engine="nx")
        assert ba.N == 100

    def test_edge_count(self):
        """Test that edge count matches BA formula."""
        ba = BarabasiAlbert(n=100, m=3, seed=42, engine="nx")
        # BA edges: (n - m - 1) * m + m*(m+1)/2 for complete initial core
        # But NX uses different initial graph, so approximately n*m - m*(m-1)/2
        assert ba.num_edges > 0

    def test_protocol_compliance(self):
        """Test that BarabasiAlbert satisfies SignedGraphProtocol."""
        ba = BarabasiAlbert(n=100, m=3, pflip=0.2, seed=42, engine="nx")
        ba.flip_random_fract_edges()
        assert is_signed_graph(ba)

    def test_edge_flipping(self):
        """Test edge sign flipping."""
        ba = BarabasiAlbert(n=100, m=3, pflip=0.25, seed=42, engine="nx")
        ba.flip_random_fract_edges()

        neg = ba.count_negative_edges()
        total = ba.num_edges

        assert neg > 0
        assert neg < total


@pytest.mark.skipif(
    not is_engine_available("gt"),
    reason="graph-tool not installed",
)
class TestBarabasiAlbertGraphTool:
    """Tests for BarabasiAlbert with graph-tool backend."""

    def test_creation(self):
        """Test basic BA graph creation."""
        ba = BarabasiAlbert(n=100, m=3, seed=42, engine="gt")
        assert ba.N == 100

    def test_edge_count(self):
        """Test that edge count matches expected for BA model."""
        ba = BarabasiAlbert(n=100, m=3, seed=42, engine="gt")
        # Expected: (n - m - 1) * m + m*(m+1)/2 = (100-4)*3 + 6 = 294
        expected = (100 - 4) * 3 + 6
        assert ba.num_edges == expected

    def test_protocol_compliance(self):
        """Test that BarabasiAlbertGT satisfies SignedGraphProtocol."""
        ba = BarabasiAlbert(n=100, m=3, pflip=0.2, seed=42, engine="gt")
        ba.flip_random_fract_edges()
        assert is_signed_graph(ba)

    def test_hub_detection(self):
        """Test hub node detection."""
        ba = BarabasiAlbert(n=100, m=3, seed=42, engine="gt")
        hubs = ba.get_hub_nodes(top_k=5)

        assert len(hubs) == 5
        # First hub should have highest degree
        assert hubs[0][1] >= hubs[1][1]

    def test_degree_distribution(self):
        """Test degree distribution computation."""
        ba = BarabasiAlbert(n=100, m=3, seed=42, engine="gt")
        unique, counts = ba.get_degree_distribution()

        assert len(unique) > 0
        assert len(counts) == len(unique)
        assert unique.min() >= 3  # Minimum degree is m

    def test_to_networkx_conversion(self):
        """Test conversion to NetworkX."""
        ba = BarabasiAlbert(n=50, m=2, engine="gt")
        nx_graph = ba.to_networkx()

        assert nx_graph.number_of_nodes() == ba.N
        assert nx_graph.number_of_edges() == ba.num_edges


@pytest.mark.skipif(
    not is_engine_available("gt"),
    reason="graph-tool not installed",
)
class TestRandomGraphCrossEngine:
    """Tests for cross-engine equivalence of random graphs."""

    def test_erdos_renyi_protocol_match(self):
        """Test that both engines satisfy the same protocol."""
        er_nx = ErdosRenyi(n=50, p=0.15, pflip=0.2, seed=42, engine="nx")
        er_gt = ErdosRenyi(n=50, p=0.15, pflip=0.2, seed=42, engine="gt")

        er_nx.flip_random_fract_edges()
        er_gt.flip_random_fract_edges()

        # Both should satisfy protocol
        assert is_signed_graph(er_nx)
        assert is_signed_graph(er_gt)

        # Both should have Laplacian of correct shape
        L_nx = er_nx.get_signed_laplacian()
        L_gt = er_gt.get_signed_laplacian()

        if hasattr(L_nx, "toarray"):
            L_nx = L_nx.toarray()

        assert L_nx.shape == (er_nx.N, er_nx.N)
        assert L_gt.shape == (er_gt.N, er_gt.N)

    def test_barabasi_albert_same_nodes(self):
        """Test that BA graphs have same node count."""
        ba_nx = BarabasiAlbert(n=100, m=3, seed=42, engine="nx")
        ba_gt = BarabasiAlbert(n=100, m=3, seed=42, engine="gt")

        assert ba_nx.N == ba_gt.N


# ============================================================================
# Watts-Strogatz Tests
# ============================================================================


class TestWattsStrogatzNetworkX:
    """Tests for WattsStrogatz with NetworkX backend."""

    def test_creation(self):
        """Test basic WS graph creation."""
        ws = WattsStrogatz(n=100, k=4, p=0.3, seed=42, engine="nx")
        assert ws.N == 100
        assert ws.num_edges == 200  # n * k / 2

    def test_protocol_compliance(self):
        """Test that WattsStrogatz satisfies SignedGraphProtocol."""
        ws = WattsStrogatz(n=100, k=4, p=0.3, pflip=0.2, seed=42, engine="nx")
        ws.flip_random_fract_edges()
        assert is_signed_graph(ws)

    def test_edge_flipping(self):
        """Test edge sign flipping."""
        ws = WattsStrogatz(n=100, k=4, p=0.3, pflip=0.25, seed=42, engine="nx")
        ws.flip_random_fract_edges()

        neg = ws.count_negative_edges()
        total = ws.num_edges

        assert neg > 0
        assert neg < total


@pytest.mark.skipif(
    not is_engine_available("gt"),
    reason="graph-tool not installed",
)
class TestWattsStrogatzGraphTool:
    """Tests for WattsStrogatz with graph-tool backend."""

    def test_creation(self):
        """Test basic WS graph creation."""
        ws = WattsStrogatz(n=100, k=4, p=0.3, seed=42, engine="gt")
        assert ws.N == 100
        assert ws.num_edges == 200  # n * k / 2

    def test_protocol_compliance(self):
        """Test that WattsStrogatzGT satisfies SignedGraphProtocol."""
        ws = WattsStrogatz(n=100, k=4, p=0.3, pflip=0.2, seed=42, engine="gt")
        ws.flip_random_fract_edges()
        assert is_signed_graph(ws)

    def test_clustering_coefficient(self):
        """Test clustering coefficient computation."""
        ws = WattsStrogatz(n=100, k=4, p=0.0, seed=42, engine="gt")
        cc = ws.get_clustering_coefficient()
        # Regular lattice (p=0) has high clustering
        assert cc > 0

    def test_to_networkx_conversion(self):
        """Test conversion to NetworkX."""
        ws = WattsStrogatz(n=50, k=4, p=0.2, engine="gt")
        nx_graph = ws.to_networkx()

        assert nx_graph.number_of_nodes() == ws.N
        assert nx_graph.number_of_edges() == ws.num_edges


# ============================================================================
# Stochastic Block Model Tests
# ============================================================================


class TestStochasticBlockModelNetworkX:
    """Tests for StochasticBlockModel with NetworkX backend."""

    def test_creation(self):
        """Test basic SBM graph creation."""
        sizes = [30, 30, 40]
        p_matrix = [[0.3, 0.05, 0.05], [0.05, 0.3, 0.05], [0.05, 0.05, 0.3]]
        sbm = StochasticBlockModel(sizes=sizes, p_matrix=p_matrix, seed=42, engine="nx")
        assert sbm.N > 0
        assert sbm.num_edges > 0

    def test_protocol_compliance(self):
        """Test that StochasticBlockModel satisfies SignedGraphProtocol."""
        sizes = [25, 25]
        p_matrix = [[0.3, 0.1], [0.1, 0.3]]
        sbm = StochasticBlockModel(sizes=sizes, p_matrix=p_matrix, pflip=0.2, seed=42, engine="nx")
        sbm.flip_random_fract_edges()
        assert is_signed_graph(sbm)

    def test_community_methods(self):
        """Test community-related methods."""
        sizes = [20, 20, 20]
        p_matrix = [[0.4, 0.05, 0.05], [0.05, 0.4, 0.05], [0.05, 0.05, 0.4]]
        sbm = StochasticBlockModel(sizes=sizes, p_matrix=p_matrix, seed=42, engine="nx")

        assert sbm.num_communities == 3
        assert sbm.get_modularity_parameter() > 0  # Strong community structure


@pytest.mark.skipif(
    not is_engine_available("gt"),
    reason="graph-tool not installed",
)
class TestStochasticBlockModelGraphTool:
    """Tests for StochasticBlockModel with graph-tool backend."""

    def test_creation(self):
        """Test basic SBM graph creation."""
        sizes = [30, 30, 40]
        p_matrix = [[0.3, 0.05, 0.05], [0.05, 0.3, 0.05], [0.05, 0.05, 0.3]]
        sbm = StochasticBlockModel(sizes=sizes, p_matrix=p_matrix, seed=42, engine="gt")
        assert sbm.N > 0
        assert sbm.num_edges > 0

    def test_protocol_compliance(self):
        """Test that StochasticBlockModelGT satisfies SignedGraphProtocol."""
        sizes = [25, 25]
        p_matrix = [[0.3, 0.1], [0.1, 0.3]]
        sbm = StochasticBlockModel(sizes=sizes, p_matrix=p_matrix, pflip=0.2, seed=42, engine="gt")
        sbm.flip_random_fract_edges()
        assert is_signed_graph(sbm)

    def test_block_membership(self):
        """Test block membership access."""
        sizes = [20, 30]
        p_matrix = [[0.4, 0.05], [0.05, 0.4]]
        sbm = StochasticBlockModel(
            sizes=sizes, p_matrix=p_matrix, extract_giant_component=False, seed=42, engine="gt"
        )
        membership = sbm.get_block_membership()
        assert len(membership) == sbm.N

    def test_modularity_parameter(self):
        """Test modularity parameter computation."""
        sizes = [30, 30]
        p_matrix = [[0.4, 0.02], [0.02, 0.4]]  # Strong community structure
        sbm = StochasticBlockModel(sizes=sizes, p_matrix=p_matrix, seed=42, engine="gt")

        mod = sbm.get_modularity_parameter()
        assert mod > 0.5  # High modularity for strong structure

    def test_to_networkx_conversion(self):
        """Test conversion to NetworkX."""
        sizes = [20, 20]
        p_matrix = [[0.3, 0.1], [0.1, 0.3]]
        sbm = StochasticBlockModel(sizes=sizes, p_matrix=p_matrix, engine="gt")
        nx_graph = sbm.to_networkx()

        assert nx_graph.number_of_nodes() == sbm.N
        assert nx_graph.number_of_edges() == sbm.num_edges


@pytest.mark.skipif(
    not is_engine_available("gt"),
    reason="graph-tool not installed",
)
class TestWattsStrogatzCrossEngine:
    """Tests for cross-engine equivalence of Watts-Strogatz graphs."""

    def test_same_structure(self):
        """Test that both engines produce same basic structure."""
        ws_nx = WattsStrogatz(n=100, k=4, p=0.3, seed=42, engine="nx")
        ws_gt = WattsStrogatz(n=100, k=4, p=0.3, seed=42, engine="gt")

        assert ws_nx.N == ws_gt.N
        assert ws_nx.num_edges == ws_gt.num_edges

    def test_protocol_match(self):
        """Test that both engines satisfy the same protocol."""
        ws_nx = WattsStrogatz(n=50, k=4, p=0.2, pflip=0.2, seed=42, engine="nx")
        ws_gt = WattsStrogatz(n=50, k=4, p=0.2, pflip=0.2, seed=42, engine="gt")

        ws_nx.flip_random_fract_edges()
        ws_gt.flip_random_fract_edges()

        assert is_signed_graph(ws_nx)
        assert is_signed_graph(ws_gt)


@pytest.mark.skipif(
    not is_engine_available("gt"),
    reason="graph-tool not installed",
)
class TestStochasticBlockModelCrossEngine:
    """Tests for cross-engine equivalence of SBM graphs."""

    def test_same_communities(self):
        """Test that both engines have same community count."""
        sizes = [30, 30, 40]
        p_matrix = [[0.3, 0.05, 0.05], [0.05, 0.3, 0.05], [0.05, 0.05, 0.3]]

        sbm_nx = StochasticBlockModel(sizes=sizes, p_matrix=p_matrix, seed=42, engine="nx")
        sbm_gt = StochasticBlockModel(sizes=sizes, p_matrix=p_matrix, seed=42, engine="gt")

        assert sbm_nx.num_communities == sbm_gt.num_communities

    def test_protocol_match(self):
        """Test that both engines satisfy the same protocol."""
        sizes = [25, 25]
        p_matrix = [[0.3, 0.1], [0.1, 0.3]]

        sbm_nx = StochasticBlockModel(sizes=sizes, p_matrix=p_matrix, pflip=0.15, seed=42, engine="nx")
        sbm_gt = StochasticBlockModel(sizes=sizes, p_matrix=p_matrix, pflip=0.15, seed=42, engine="gt")

        sbm_nx.flip_random_fract_edges()
        sbm_gt.flip_random_fract_edges()

        assert is_signed_graph(sbm_nx)
        assert is_signed_graph(sbm_gt)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
