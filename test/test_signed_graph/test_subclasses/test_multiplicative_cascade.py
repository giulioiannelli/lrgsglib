"""
Unit tests for MultiplicativeCascadeGraph class.
"""
import pytest
import networkx as nx
import numpy as np

from lrgsglib.graphs.nx import MultiplicativeCascadeGraphNX as MultiplicativeCascadeGraph
from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph


class TestMultiplicativeCascadeGraphCreation:
    """Test graph creation and initialization."""

    def test_default_parameters(self):
        """Test creation with default parameters."""
        msg = MultiplicativeCascadeGraph(only_const_mode=True)
        assert isinstance(msg, MultiplicativeCascadeGraph)
        assert isinstance(msg, SignedGraph)
        assert msg.only_const_mode is True

    def test_custom_parameters(self):
        """Test creation with custom parameters."""
        msg = MultiplicativeCascadeGraph(
            p1=0.8,
            p2=0.6,
            p3=0.6,
            p4=0.8,
            fraction=0.4,
            iterations=3,
            only_const_mode=True,
        )
        assert msg.p1 == 0.8
        assert msg.p2 == 0.6
        assert msg.p3 == 0.6
        assert msg.p4 == 0.8
        assert msg.fraction == 0.4
        assert msg.iterations == 3

    def test_seed_probabilities_attribute(self):
        """Test that seed_probabilities attribute is set correctly."""
        msg = MultiplicativeCascadeGraph(
            p1=0.9, p2=0.7, p3=0.5, p4=0.3, only_const_mode=True
        )
        assert msg.seed_probabilities == (0.9, 0.7, 0.5, 0.3)

    def test_sample_fraction_attribute(self):
        """Test that sample_fraction equals fraction."""
        msg = MultiplicativeCascadeGraph(fraction=0.25, only_const_mode=True)
        assert msg.sample_fraction == 0.25

    def test_variant_standard(self):
        """Test standard variant."""
        msg = MultiplicativeCascadeGraph(variant="standard", only_const_mode=True)
        assert msg.variant == "standard"

    def test_variant_exp_clocks(self):
        """Test exp_clocks variant."""
        msg = MultiplicativeCascadeGraph(variant="exp_clocks", only_const_mode=True)
        assert msg.variant == "exp_clocks"


class TestMultiplicativeCascadeGraphGeneration:
    """Test graph generation."""

    def test_graph_generation_standard(self):
        """Test that graph is generated with standard variant."""
        msg = MultiplicativeCascadeGraph(
            p1=0.8, p2=0.6, p3=0.6, p4=0.8, fraction=0.5, iterations=2
        )
        assert isinstance(msg.G, nx.Graph)
        assert msg.N > 0
        assert msg.probability_matrix is not None

    def test_graph_generation_exp_clocks(self):
        """Test that graph is generated with exp_clocks variant."""
        msg = MultiplicativeCascadeGraph(
            p1=0.8,
            p2=0.6,
            p3=0.6,
            p4=0.8,
            fraction=0.5,
            iterations=2,
            variant="exp_clocks",
        )
        assert isinstance(msg.G, nx.Graph)
        assert msg.N > 0
        assert msg.probability_matrix is not None

    def test_graph_not_generated_in_const_mode(self):
        """Test that graph is not generated in only_const_mode."""
        msg = MultiplicativeCascadeGraph(only_const_mode=True)
        assert msg.G.number_of_nodes() == 0
        assert msg.probability_matrix is None

    def test_probability_matrix_shape(self):
        """Test that probability matrix has correct shape."""
        iterations = 3
        msg = MultiplicativeCascadeGraph(iterations=iterations)
        # Size is 2^(iterations+1) based on implementation
        expected_size = 2**(iterations + 1)
        assert msg.probability_matrix.shape == (expected_size, expected_size)

    def test_periodic_boundary_conditions(self):
        """Test graph generation with periodic boundaries."""
        msg = MultiplicativeCascadeGraph(
            p1=0.8, p2=0.6, p3=0.6, p4=0.8, fraction=0.5, iterations=2, periodic=True
        )
        assert isinstance(msg.G, nx.Graph)
        assert msg.periodic is True

    def test_stochastic_mode(self):
        """Test graph generation in stochastic mode."""
        msg = MultiplicativeCascadeGraph(
            p1=0.8,
            p2=0.6,
            p3=0.6,
            p4=0.8,
            fraction=0.5,
            iterations=2,
            stochastic=True,
        )
        assert isinstance(msg.G, nx.Graph)
        assert msg.stochastic is True


class TestMultiplicativeCascadeGraphNodeCount:
    """Test expected node count calculations."""

    def test_expected_num_nodes_const_mode(self):
        """Test that expected_num_nodes returns 0 in const mode."""
        msg = MultiplicativeCascadeGraph(only_const_mode=True)
        assert msg.get_expected_num_nodes() == 0

    def test_expected_num_nodes_formula(self):
        """Test expected node count calculation."""
        iterations = 3
        fraction = 0.5
        msg = MultiplicativeCascadeGraph(iterations=iterations, fraction=fraction)
        # Size is 2^(iterations+1) based on implementation
        size = 2**(iterations + 1)
        expected = max(1, int(round(fraction * (size**2))))
        assert msg.get_expected_num_nodes() == expected

    def test_actual_vs_expected_nodes(self):
        """Test that actual node count is reasonable."""
        msg = MultiplicativeCascadeGraph(
            p1=0.8, p2=0.6, p3=0.6, p4=0.8, fraction=0.5, iterations=3
        )
        expected = msg.get_expected_num_nodes()
        actual = msg.N
        # Actual N can vary significantly due to probabilistic generation
        # Just check that we got some nodes
        assert actual > 0
        assert actual <= expected * 2  # Reasonable upper bound


class TestMultiplicativeCascadeGraphRepresentation:
    """Test string representations."""

    def test_repr(self):
        """Test __repr__ method."""
        msg = MultiplicativeCascadeGraph(
            p1=0.8, p2=0.6, p3=0.6, p4=0.8, fraction=0.5, iterations=2
        )
        repr_str = repr(msg)
        assert "MultiplicativeCascadeGraph" in repr_str
        assert "N=" in repr_str
        assert "variant=" in repr_str

    def test_repr_variant_included(self):
        """Test that variant is included in repr."""
        msg = MultiplicativeCascadeGraph(
            variant="exp_clocks", iterations=2, only_const_mode=True
        )
        repr_str = repr(msg)
        assert "exp_clocks" in repr_str


class TestMultiplicativeCascadeGraphIntegration:
    """Test integration with SignedGraph methods."""

    def test_inherits_from_signed_graph(self):
        """Test that MultiplicativeCascadeGraph inherits from SignedGraph."""
        msg = MultiplicativeCascadeGraph(only_const_mode=True)
        assert isinstance(msg, SignedGraph)

    def test_laplacian_matrix_computation(self):
        """Test that Laplacian matrix can be computed."""
        msg = MultiplicativeCascadeGraph(
            p1=0.8, p2=0.6, p3=0.6, p4=0.8, fraction=0.5, iterations=2
        )
        lap = msg.lap
        assert lap is not None
        assert lap.shape[0] == msg.N
        assert lap.shape[1] == msg.N

    def test_spectral_computation(self):
        """Test that spectral methods work."""
        msg = MultiplicativeCascadeGraph(
            p1=0.8, p2=0.6, p3=0.6, p4=0.8, fraction=0.5, iterations=2
        )
        msg.compute_laplacian_spectrum()
        assert msg.eigv is not None
        assert len(msg.eigv) == msg.N

    def test_properties_N_Ne(self):
        """Test that N and Ne properties work."""
        msg = MultiplicativeCascadeGraph(
            p1=0.8, p2=0.6, p3=0.6, p4=0.8, fraction=0.5, iterations=2
        )
        assert isinstance(msg.N, int)
        assert msg.N > 0
        assert isinstance(msg.Ne, int)
        assert msg.Ne >= 0

    def test_pflip_parameter(self):
        """Test that pflip parameter is respected."""
        msg = MultiplicativeCascadeGraph(
            p1=0.8, p2=0.6, p3=0.6, p4=0.8, fraction=0.5, iterations=2, pflip=0.5
        )
        assert msg.pflip == 0.5

    def test_sgpathn_attribute(self):
        """Test that sgpathn is set correctly."""
        msg = MultiplicativeCascadeGraph(only_const_mode=True)
        assert msg.sgpathn == "msg_mc"

    def test_custom_sgpathn(self):
        """Test that custom sgpathn can be set."""
        custom_path = "custom_mc_path"
        msg = MultiplicativeCascadeGraph(sgpathn=custom_path, only_const_mode=True)
        assert msg.sgpathn == custom_path


class TestMultiplicativeCascadeGraphEdgeCases:
    """Test edge cases and error conditions."""

    def test_small_fraction(self):
        """Test with small fraction."""
        msg = MultiplicativeCascadeGraph(fraction=0.1, iterations=2)
        # Should create at least some nodes
        assert msg.N >= 1

    def test_full_fraction(self):
        """Test with fraction=1.0 (full graph)."""
        iterations = 2
        msg = MultiplicativeCascadeGraph(fraction=1.0, iterations=iterations)
        # Size is 2^(iterations+1)
        size = 2**(iterations + 1)
        expected = size**2
        # Should be close to full size
        assert msg.N > 0
        assert msg.N <= expected * 1.5  # Allow some margin

    def test_single_iteration(self):
        """Test with single iteration."""
        msg = MultiplicativeCascadeGraph(iterations=1, fraction=0.5)
        assert isinstance(msg.G, nx.Graph)
        assert msg.N > 0

    def test_multiple_iterations(self):
        """Test with multiple iterations."""
        msg = MultiplicativeCascadeGraph(iterations=4, fraction=0.3)
        assert isinstance(msg.G, nx.Graph)
        assert msg.N > 0
