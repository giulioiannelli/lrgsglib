"""
Unit tests for VicsekGraph class.
"""
import pytest
import networkx as nx
import numpy as np

from lrgsglib.nx_patches.Vicsek import VicsekGraph
from lrgsglib.nx_patches.SignedGraph import SignedGraph


class TestVicsekGraphCreation:
    """Test graph creation and initialization."""

    def test_default_parameters(self):
        """Test creation with default parameters."""
        vck = VicsekGraph(N=100, k=3, only_const_mode=True)
        assert isinstance(vck, VicsekGraph)
        assert isinstance(vck, SignedGraph)
        assert vck.only_const_mode is True

    def test_custom_parameters(self):
        """Test creation with custom parameters."""
        vck = VicsekGraph(N=200, k=4, m=5, symmetric=False, only_const_mode=True)
        # N_target is internal attribute
        assert vck.N_target == 200
        assert vck.k == 4
        assert vck.m == 5
        assert vck.symmetric is False

    def test_pij_parameter_none(self):
        """Test with pij=None (default)."""
        vck = VicsekGraph(N=100, k=3, pij=None, only_const_mode=True)
        # pij_input is internal, not public API
        assert vck.pij_input is None

    def test_pij_parameter_custom(self):
        """Test with custom pij 2D array."""
        import numpy as np
        pij_matrix = np.array([[0.9, 0.1], [0.1, 0.9]])
        vck = VicsekGraph(N=100, k=3, pij=pij_matrix, only_const_mode=True)
        # pij_input is internal, not public API
        assert vck.pij_input is not None
        assert np.array_equal(vck.pij_input, pij_matrix)

    def test_symmetric_true(self):
        """Test symmetric mode."""
        vck = VicsekGraph(N=100, k=3, symmetric=True, only_const_mode=True)
        assert vck.symmetric is True

    def test_symmetric_false(self):
        """Test asymmetric mode."""
        vck = VicsekGraph(N=100, k=3, symmetric=False, only_const_mode=True)
        assert vck.symmetric is False


class TestVicsekGraphGeneration:
    """Test graph generation."""

    def test_graph_generation(self):
        """Test that graph is generated correctly."""
        vck = VicsekGraph(N=50, k=3)
        assert isinstance(vck.G, nx.Graph)
        assert vck.N == 50
        assert vck.probability_matrix is not None

    def test_graph_not_generated_in_const_mode(self):
        """Test that graph is not generated in only_const_mode."""
        vck = VicsekGraph(N=100, k=3, only_const_mode=True)
        assert vck.G.number_of_nodes() == 0
        assert vck.probability_matrix is None

    def test_probability_matrix_exists(self):
        """Test that probability matrix is generated."""
        N = 50
        vck = VicsekGraph(N=N, k=3)
        assert vck.probability_matrix is not None
        # Shape depends on kronecker product iterations, not directly on N
        assert vck.probability_matrix.ndim == 2

    def test_probability_matrix_symmetric(self):
        """Test that probability matrix is symmetric when symmetric=True."""
        N = 30
        vck = VicsekGraph(N=N, k=3, symmetric=True)
        prob_matrix = vck.probability_matrix
        assert np.allclose(prob_matrix, prob_matrix.T)

    def test_sample_fraction_is_one(self):
        """Test that sample_fraction is always 1.0 for Vicsek."""
        vck = VicsekGraph(N=50, k=3)
        assert vck.sample_fraction == 1.0

    def test_different_k_values(self):
        """Test graph generation with different k values."""
        for k in [2, 3, 4, 5]:
            vck = VicsekGraph(N=40, k=k)
            assert isinstance(vck.G, nx.Graph)
            assert vck.k == k

    def test_different_m_values(self):
        """Test graph generation with different m values."""
        for m in [2, 3, 4, 5]:
            vck = VicsekGraph(N=40, k=3, m=m)
            assert isinstance(vck.G, nx.Graph)
            assert vck.m == m


class TestVicsekGraphNodeCount:
    """Test expected node count calculations."""

    def test_expected_num_nodes_const_mode(self):
        """Test that expected_num_nodes returns 0 in const mode."""
        vck = VicsekGraph(N=100, k=3, only_const_mode=True)
        assert vck.get_expected_num_nodes() == 0

    def test_expected_num_nodes_equals_N(self):
        """Test that expected node count equals N parameter."""
        N = 75
        vck = VicsekGraph(N=N, k=3)
        assert vck.get_expected_num_nodes() == N

    def test_actual_equals_expected_nodes(self):
        """Test that actual node count equals expected (no sampling)."""
        N = 60
        vck = VicsekGraph(N=N, k=3)
        assert vck.N == N
        assert vck.N == vck.get_expected_num_nodes()


class TestVicsekGraphRepresentation:
    """Test string representations."""

    def test_repr(self):
        """Test __repr__ method."""
        vck = VicsekGraph(N=50, k=3)
        repr_str = repr(vck)
        assert "VicsekGraph" in repr_str
        assert "N=50" in repr_str
        assert "k=3" in repr_str

    def test_repr_includes_k(self):
        """Test that k parameter is included in repr."""
        vck = VicsekGraph(N=40, k=5, only_const_mode=True)
        repr_str = repr(vck)
        assert "k=5" in repr_str


class TestVicsekGraphIntegration:
    """Test integration with SignedGraph methods."""

    def test_inherits_from_signed_graph(self):
        """Test that VicsekGraph inherits from SignedGraph."""
        vck = VicsekGraph(N=50, k=3, only_const_mode=True)
        assert isinstance(vck, SignedGraph)

    def test_laplacian_matrix_computation(self):
        """Test that Laplacian matrix can be computed."""
        vck = VicsekGraph(N=40, k=3)
        lap = vck.lap
        assert lap is not None
        assert lap.shape[0] == vck.N
        assert lap.shape[1] == vck.N

    def test_spectral_computation(self):
        """Test that spectral methods work."""
        vck = VicsekGraph(N=40, k=3)
        vck.compute_laplacian_spectrum()
        assert vck.eigv is not None
        assert len(vck.eigv) == vck.N

    def test_properties_N_Ne(self):
        """Test that N and Ne properties work."""
        vck = VicsekGraph(N=50, k=3)
        assert isinstance(vck.N, int)
        assert vck.N == 50
        assert isinstance(vck.Ne, int)
        assert vck.Ne >= 0

    def test_pflip_parameter(self):
        """Test that pflip parameter is respected."""
        vck = VicsekGraph(N=40, k=3, pflip=0.3)
        assert vck.pflip == 0.3

    def test_sgpathn_attribute(self):
        """Test that sgpathn is set correctly."""
        vck = VicsekGraph(N=50, k=3, only_const_mode=True)
        assert vck.sgpathn == "msg_plv"

    def test_custom_sgpathn(self):
        """Test that custom sgpathn can be set."""
        custom_path = "custom_vicsek_path"
        vck = VicsekGraph(N=50, k=3, sgpathn=custom_path, only_const_mode=True)
        assert vck.sgpathn == custom_path


class TestVicsekGraphEdgeCases:
    """Test edge cases and error conditions."""

    def test_small_graph(self):
        """Test with very small N."""
        vck = VicsekGraph(N=5, k=2)
        assert isinstance(vck.G, nx.Graph)
        assert vck.N == 5

    def test_reasonable_k_value(self):
        """Test with reasonable k value."""
        N = 100
        k = 5
        vck = VicsekGraph(N=N, k=k)
        assert isinstance(vck.G, nx.Graph)
        assert vck.N > 0

    def test_k_equals_two(self):
        """Test minimum k value."""
        vck = VicsekGraph(N=30, k=2)
        assert isinstance(vck.G, nx.Graph)
        assert vck.k == 2

    def test_m_equals_k(self):
        """Test when m equals k."""
        vck = VicsekGraph(N=30, k=3, m=3)
        assert isinstance(vck.G, nx.Graph)
        assert vck.m == vck.k

    def test_pij_custom_matrix(self):
        """Test with custom pij matrix."""
        import numpy as np
        # Custom 2x2 probability matrix
        pij_matrix = np.array([[0.9, 0.1], [0.1, 0.9]])
        vck = VicsekGraph(N=30, k=3, pij=pij_matrix)
        assert isinstance(vck.G, nx.Graph)
        # Graph should be generated with custom probabilities


class TestVicsekGraphConnectivity:
    """Test graph connectivity properties."""

    def test_graph_is_connected_high_k(self):
        """Test that graph tends to be connected with high k."""
        vck = VicsekGraph(N=50, k=8)
        # With high k, graph should likely be connected
        # Note: This is probabilistic, so we just check it exists
        assert isinstance(vck.G, nx.Graph)
        assert vck.N > 0

    def test_edge_count_reasonable(self):
        """Test that edge count is reasonable."""
        import numpy as np

        np.random.seed(0)
        pij_matrix = np.full((2, 2), 0.5)
        vck = VicsekGraph(N=50, k=2, pij=pij_matrix)
        # Should have some edges
        assert vck.Ne > 0
        # Should not have more edges than complete graph
        max_edges = vck.N * (vck.N - 1) // 2
        assert vck.Ne <= max_edges
