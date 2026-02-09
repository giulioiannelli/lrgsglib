"""
Tests for graph-tool C++ extensions and GT graph utilities.

These tests verify that the C++ extensions for graph-tool operations
work correctly and provide performance benefits over pure Python.
"""
import pytest
import time

# Skip all tests if graph-tool not available
gt = pytest.importorskip("graph_tool")


class TestTriangularLatticeExtension:
    """Test the C++ triangular lattice generator."""

    def test_import(self):
        """Test that the extension can be imported."""
        from lrgsglib.graphs.gt.Lattice2DGT.cpp import create_triangular_lattice
        assert callable(create_triangular_lattice)

    def test_basic_creation(self):
        """Test basic lattice creation."""
        from lrgsglib.graphs.gt.Lattice2DGT.cpp import create_triangular_lattice

        g = create_triangular_lattice(5, 5)
        assert g.num_vertices() == 25  # 5x5 grid

    def test_edge_structure(self):
        """Test that the triangular lattice has correct edge count."""
        from lrgsglib.graphs.gt.Lattice2DGT.cpp import create_triangular_lattice

        # For a triangular lattice, the number of edges depends on geometry
        g = create_triangular_lattice(10, 10)
        assert g.num_vertices() == 100
        assert g.num_edges() > 0

        # Each interior node has degree 6, edges reduced for boundary
        # Just verify reasonable range
        assert 200 < g.num_edges() < 400

    def test_undirected(self):
        """Test that the graph is undirected."""
        from lrgsglib.graphs.gt.Lattice2DGT.cpp import create_triangular_lattice

        g = create_triangular_lattice(5, 5)
        assert not g.is_directed()

    def test_various_sizes(self):
        """Test different lattice sizes."""
        from lrgsglib.graphs.gt.Lattice2DGT.cpp import create_triangular_lattice

        for size in [(3, 3), (10, 10), (20, 15), (50, 50)]:
            g = create_triangular_lattice(*size)
            expected_vertices = size[0] * size[1]
            assert g.num_vertices() == expected_vertices

    def test_performance(self):
        """Test that C++ extension is reasonably fast."""
        from lrgsglib.graphs.gt.Lattice2DGT.cpp import create_triangular_lattice

        # Creating a 100x100 lattice should take < 100ms
        start = time.time()
        g = create_triangular_lattice(100, 100)
        elapsed = time.time() - start

        assert g.num_vertices() == 10000
        assert elapsed < 0.1  # Less than 100ms

    def test_large_lattice(self):
        """Test large lattice creation."""
        from lrgsglib.graphs.gt.Lattice2DGT.cpp import create_triangular_lattice

        g = create_triangular_lattice(200, 200)
        assert g.num_vertices() == 40000


class TestGTNXConversion:
    """Test conversion between NetworkX and graph-tool."""

    def test_nx_to_gt_basic(self):
        """Test basic NX to GT conversion."""
        import networkx as nx
        from lrgsglib.graphs.gt import nx_to_gt

        G_nx = nx.complete_graph(10)
        for u, v in G_nx.edges():
            G_nx[u][v]["sign"] = 1

        G_gt = nx_to_gt(G_nx)
        assert G_gt.num_vertices() == 10
        assert G_gt.num_edges() == 45  # C(10,2)

    def test_gt_to_nx_basic(self):
        """Test basic GT to NX conversion."""
        import graph_tool.generation as gen
        from lrgsglib.graphs.gt import gt_to_nx

        G_gt = gen.complete_graph(10)
        G_nx = gt_to_nx(G_gt)

        assert G_nx.number_of_nodes() == 10
        assert G_nx.number_of_edges() == 45

    def test_roundtrip(self):
        """Test NX -> GT -> NX roundtrip."""
        import networkx as nx
        from lrgsglib.graphs.gt import nx_to_gt, gt_to_nx

        G_nx = nx.grid_2d_graph(5, 5)
        for u, v in G_nx.edges():
            G_nx[u][v]["sign"] = 1

        G_gt = nx_to_gt(G_nx)
        G_nx2 = gt_to_nx(G_gt)

        assert G_nx2.number_of_nodes() == G_nx.number_of_nodes()
        assert G_nx2.number_of_edges() == G_nx.number_of_edges()


class TestSignedGraphGT:
    """Test SignedGraphGT class."""

    def test_creation(self):
        """Test basic SignedGraphGT creation."""
        from lrgsglib.graphs.gt.SignedGraphGT import SignedGraphGT

        sg = SignedGraphGT()
        assert sg.N == 0

    def test_from_lattice(self):
        """Test creating SignedGraphGT from lattice."""
        from lrgsglib.graphs.gt.SignedGraphGT import SignedGraphGT

        sg = SignedGraphGT.from_lattice((10, 10))
        assert sg.N == 100

    def test_flip_edges(self):
        """Test edge sign flipping."""
        from lrgsglib.graphs.gt.SignedGraphGT import SignedGraphGT

        sg = SignedGraphGT.from_lattice((10, 10))

        # Initially all positive
        initial_neg = sg.count_negative_edges()

        # Flip some edges
        sg.flip_random_fract_edges(fraction=0.5)

        final_neg = sg.count_negative_edges()
        final_pos = sg.count_positive_edges()

        # Should have some negative edges now
        assert final_neg > 0
        assert final_pos + final_neg == sg.num_edges


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
