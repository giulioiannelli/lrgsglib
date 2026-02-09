"""
Tests for the 3 critical bugs fixed in Phase 1.

Bug 1: Missing nodes_in() method
Bug 2: Broken unflip_all() method
Bug 3: Inconsistent attribute naming (specific_heat vs entropy_derivative)
"""
import pytest
import networkx as nx
import warnings


class TestBug1NodesIn:
    """Tests for Bug 1: Missing nodes_in() method."""

    def test_nodes_in_exists(self, small_graph):
        """Test that nodes_in() method exists on SignedGraph."""
        from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph

        sg = SignedGraph(G=small_graph)
        assert hasattr(sg, 'nodes_in'), "nodes_in() method should exist"

    def test_nodes_in_returns_correct_nodes(self, small_graph):
        """Test that nodes_in() returns the correct list of nodes."""
        from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph

        sg = SignedGraph(G=small_graph)
        nodes = sg.nodes_in()

        assert isinstance(nodes, list), "nodes_in() should return a list"
        assert len(nodes) == small_graph.number_of_nodes()
        assert set(nodes) == set(small_graph.nodes())

    def test_nodes_in_with_multiple_representations(self, small_graph):
        """Test nodes_in() with different graph representations."""
        from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph
        from lrgsglib.config.const import SG_REPR

        sg = SignedGraph(G=small_graph)
        nodes = sg.nodes_in(on_g=SG_REPR)

        assert len(nodes) == small_graph.number_of_nodes()


class TestBug2UnflipAll:
    """Tests for Bug 2: Broken unflip_all() method."""

    def test_unflip_all_basic(self, small_graph, random_seed):
        """Test that unflip_all() runs without error."""
        from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph

        sg = SignedGraph(G=small_graph, pflip=0.5, seed=random_seed)

        # Should not raise an error
        sg.unflip_all()

    def test_unflip_all_restores_positive(self, small_graph, random_seed):
        """Test that unflip_all() flips all negative edges back to positive."""
        from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph

        sg = SignedGraph(G=small_graph, pflip=0.5, seed=random_seed)

        # Before unflip_all, should have negative edges
        neg_edges_before = len(sg.fleset[sg.on_g])
        assert neg_edges_before > 0, "Should have negative edges with pflip=0.5"

        # After unflip_all, should have no negative edges
        sg.unflip_all()
        neg_edges_after = len(sg.fleset[sg.on_g])

        assert neg_edges_after == 0, "All edges should be positive after unflip_all()"

    def test_unflip_all_updates_fleset(self, small_graph, random_seed):
        """Test that unflip_all() properly updates the fleset."""
        from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph

        sg = SignedGraph(G=small_graph, pflip=0.3, seed=random_seed)

        sg.unflip_all()

        # fleset should be empty
        assert len(sg.fleset[sg.on_g]) == 0

        # All edges should have weight +1
        for u, v, data in sg.G.edges(data=True):
            assert data['weight'] == 1, f"Edge ({u}, {v}) should have weight +1"


class TestBug3AttributeNaming:
    """Tests for Bug 3: Inconsistent attribute naming (specific_heat vs entropy_derivative)."""

    def test_specific_heat_attribute_exists(self, small_graph):
        """Test that specific_heat is set after entropy computation."""
        from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph

        sg = SignedGraph(G=small_graph)
        sg.compute_laplacian_spectrum_weigV()

        # Import the function
        from lrgsglib.graphs.nx.SignedGraphNX._infotheory import compute_signed_laplacian_entropy

        compute_signed_laplacian_entropy(sg)

        assert hasattr(sg, 'specific_heat'), "specific_heat attribute should exist"
        assert sg.specific_heat is not None, "specific_heat should be computed"

    def test_entropy_derivative_deprecated_alias_works(self, small_graph):
        """Test that entropy_derivative still works as deprecated alias."""
        from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph

        sg = SignedGraph(G=small_graph)
        sg.compute_laplacian_spectrum_weigV()

        from lrgsglib.graphs.nx.SignedGraphNX._infotheory import compute_signed_laplacian_entropy
        compute_signed_laplacian_entropy(sg)

        # Should have both attributes
        assert hasattr(sg, 'specific_heat')
        assert hasattr(sg, 'entropy_derivative'), "entropy_derivative deprecated alias should exist"

        # They should be the same
        import numpy as np
        np.testing.assert_array_equal(sg.specific_heat, sg.entropy_derivative)

    def test_get_specific_heat_computes_if_missing(self, small_graph):
        """Test that get_specific_heat() computes entropy if not already computed."""
        from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph

        sg = SignedGraph(G=small_graph)
        sg.compute_laplacian_spectrum_weigV()

        # get_specific_heat should compute it automatically
        result = sg.get_specific_heat()

        assert result is not None
        assert hasattr(sg, 'specific_heat')
        assert sg.specific_heat is not None

    def test_entropy_derivative_deprecation_warning(self, small_graph):
        """Test that using entropy_derivative getter shows deprecation warning."""
        from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph

        sg = SignedGraph(G=small_graph)
        sg.compute_laplacian_spectrum_weigV()

        # Using get_entropy_derivative should warn
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            result = sg.get_entropy_derivative()

            assert len(w) >= 1
            assert issubclass(w[0].category, DeprecationWarning)
            assert "specific_heat" in str(w[0].message).lower()
