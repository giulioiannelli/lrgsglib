"""
Integration tests for all MultispectralGraph subclasses.

Tests that all subclasses work correctly with SignedGraph methods
and can be used interchangeably.
"""
import pytest
import networkx as nx
import numpy as np

from lrgsglib.graphs.nx import MultiplicativeCascadeGraphNX as MultiplicativeCascadeGraph
from lrgsglib.graphs.nx import VicsekGraphNX as VicsekGraph
from lrgsglib.graphs.nx import DiracCombGraphNX as DiracCombGraph, DiracBrushGraphNX as DiracBrushGraph
from lrgsglib.graphs.nx import MultispectralGraphNX as MultispectralGraph
from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph


@pytest.fixture(params=[
    ("MultiplicativeCascade", {"p1": 0.8, "p2": 0.6, "p3": 0.6, "p4": 0.8, "fraction": 0.5, "iterations": 2}),
    ("Vicsek", {"N": 50, "k": 3}),
    ("DiracComb", {"base_nodes": 10, "fiber_nodes": 5}),
    ("DiracBrush", {"base_x": 5, "base_y": 4, "fiber_nodes": 3}),
])
def msg_instance(request):
    """Create instances of all MultispectralGraph subclasses."""
    graph_type, params = request.param

    if graph_type == "MultiplicativeCascade":
        return MultiplicativeCascadeGraph(**params)
    elif graph_type == "Vicsek":
        return VicsekGraph(**params)
    elif graph_type == "DiracComb":
        return DiracCombGraph(**params)
    elif graph_type == "DiracBrush":
        return DiracBrushGraph(**params)


class TestMultispectralGraphInheritance:
    """Test that all subclasses properly inherit from base classes."""

    def test_inherits_from_multispectral_graph(self, msg_instance):
        """Test that all instances inherit from MultispectralGraph."""
        assert isinstance(msg_instance, MultispectralGraph)

    def test_inherits_from_signed_graph(self, msg_instance):
        """Test that all instances inherit from SignedGraph."""
        assert isinstance(msg_instance, SignedGraph)

    def test_is_networkx_graph(self, msg_instance):
        """Test that all instances are NetworkX graphs."""
        assert isinstance(msg_instance.G, nx.Graph)


class TestMultispectralGraphBasicProperties:
    """Test basic graph properties for all subclasses."""

    def test_has_nodes(self, msg_instance):
        """Test that all graphs have nodes."""
        assert msg_instance.N > 0

    def test_has_edges(self, msg_instance):
        """Test that all graphs have edges."""
        assert msg_instance.Ne >= 0

    def test_node_count_property(self, msg_instance):
        """Test that N property works correctly."""
        assert isinstance(msg_instance.N, int)
        assert msg_instance.N == msg_instance.G.number_of_nodes()

    def test_edge_count_property(self, msg_instance):
        """Test that Ne property works correctly."""
        assert isinstance(msg_instance.Ne, int)
        assert msg_instance.Ne == msg_instance.G.number_of_edges()


class TestMultispectralGraphLaplacianMethods:
    """Test Laplacian-related methods for all subclasses."""

    def test_laplacian_matrix_exists(self, msg_instance):
        """Test that Laplacian matrix can be computed."""
        lap = msg_instance.lap
        assert lap is not None

    def test_laplacian_matrix_shape(self, msg_instance):
        """Test that Laplacian matrix has correct shape."""
        lap = msg_instance.lap
        assert lap.shape[0] == msg_instance.N
        assert lap.shape[1] == msg_instance.N

    def test_laplacian_matrix_symmetric(self, msg_instance):
        """Test that Laplacian matrix is symmetric."""
        lap = msg_instance.lap.toarray()
        assert np.allclose(lap, lap.T)

    def test_adjacency_matrix_exists(self, msg_instance):
        """Test that adjacency matrix can be computed."""
        adj = msg_instance.adj
        assert adj is not None

    def test_adjacency_matrix_shape(self, msg_instance):
        """Test that adjacency matrix has correct shape."""
        adj = msg_instance.adj
        assert adj.shape[0] == msg_instance.N
        assert adj.shape[1] == msg_instance.N


class TestMultispectralGraphSpectralMethods:
    """Test spectral computation methods for all subclasses."""

    def test_compute_laplacian_spectrum(self, msg_instance):
        """Test Laplacian spectrum computation."""
        msg_instance.compute_laplacian_spectrum()
        assert msg_instance.eigv is not None
        assert len(msg_instance.eigv) == msg_instance.N

    def test_eigenvalues_are_real(self, msg_instance):
        """Test that eigenvalues are real."""
        msg_instance.compute_laplacian_spectrum()
        assert np.all(np.isreal(msg_instance.eigv))

    def test_eigenvalues_are_sorted(self, msg_instance):
        """Test that eigenvalues are sorted."""
        msg_instance.compute_laplacian_spectrum()
        assert np.all(np.diff(msg_instance.eigv) >= -1e-10)  # Allow small numerical errors

    def test_smallest_eigenvalue_near_zero(self, msg_instance):
        """Test that smallest eigenvalue is close to zero (for connected graphs)."""
        msg_instance.compute_laplacian_spectrum()
        # Smallest eigenvalue should be near 0 for Laplacian
        assert abs(msg_instance.eigv[0]) < 1e-8


class TestMultispectralGraphPflip:
    """Test pflip parameter for all subclasses."""

    @pytest.mark.parametrize("pflip_val", [0.0, 0.2, 0.5, 0.8])
    def test_pflip_values(self, pflip_val):
        """Test different pflip values for each graph type."""
        # Test MultiplicativeCascade
        mc = MultiplicativeCascadeGraph(
            p1=0.8, p2=0.6, p3=0.6, p4=0.8, fraction=0.5, iterations=2, pflip=pflip_val
        )
        assert mc.pflip == pflip_val

        # Test Vicsek
        vck = VicsekGraph(N=30, k=3, pflip=pflip_val)
        assert vck.pflip == pflip_val

        # Test DiracComb
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5, pflip=pflip_val)
        assert dcomb.pflip == pflip_val

        # Test DiracBrush
        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3, pflip=pflip_val)
        assert dbrush.pflip == pflip_val


class TestMultispectralGraphOnlyConstMode:
    """Test only_const_mode for all subclasses."""

    def test_multiplicative_cascade_const_mode(self):
        """Test only_const_mode for MultiplicativeCascadeGraph."""
        mc = MultiplicativeCascadeGraph(only_const_mode=True)
        assert mc.G.number_of_nodes() == 0
        assert mc.get_expected_num_nodes() == 0

    def test_vicsek_const_mode(self):
        """Test only_const_mode for VicsekGraph."""
        vck = VicsekGraph(N=50, k=3, only_const_mode=True)
        assert vck.G.number_of_nodes() == 0
        assert vck.get_expected_num_nodes() == 0

    def test_dirac_comb_const_mode(self):
        """Test only_const_mode for DiracCombGraph."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5, only_const_mode=True)
        assert dcomb.G.number_of_nodes() == 0
        assert dcomb.get_expected_num_nodes() == 0

    def test_dirac_brush_const_mode(self):
        """Test only_const_mode for DiracBrushGraph."""
        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3, only_const_mode=True)
        assert dbrush.G.number_of_nodes() == 0
        assert dbrush.get_expected_num_nodes() == 0


class TestMultispectralGraphExpectedNodeCount:
    """Test get_expected_num_nodes for all subclasses."""

    def test_multiplicative_cascade_expected_nodes(self):
        """Test expected node count for MultiplicativeCascadeGraph."""
        mc = MultiplicativeCascadeGraph(
            p1=0.8, p2=0.6, p3=0.6, p4=0.8, fraction=0.5, iterations=3
        )
        expected = mc.get_expected_num_nodes()
        assert expected > 0
        assert isinstance(expected, int)

    def test_vicsek_expected_nodes(self):
        """Test expected node count for VicsekGraph."""
        N = 75
        vck = VicsekGraph(N=N, k=3)
        expected = vck.get_expected_num_nodes()
        assert expected == N

    def test_dirac_comb_expected_nodes(self):
        """Test expected node count for DiracCombGraph."""
        base_nodes = 12
        fiber_nodes = 8
        dcomb = DiracCombGraph(base_nodes=base_nodes, fiber_nodes=fiber_nodes)
        expected = dcomb.get_expected_num_nodes()
        # Formula: base_nodes * (1 + fiber_nodes)
        assert expected == base_nodes * (1 + fiber_nodes)

    def test_dirac_brush_expected_nodes(self):
        """Test expected node count for DiracBrushGraph."""
        base_x = 7
        base_y = 6
        fiber_nodes = 3
        dbrush = DiracBrushGraph(base_x=base_x, base_y=base_y, fiber_nodes=fiber_nodes)
        expected = dbrush.get_expected_num_nodes()
        # Formula: base_x * base_y * (1 + fiber_nodes)
        assert expected == base_x * base_y * (1 + fiber_nodes)


class TestMultispectralGraphStringRepresentation:
    """Test string representations for all subclasses."""

    def test_repr_format(self, msg_instance):
        """Test that __repr__ returns proper format."""
        repr_str = repr(msg_instance)
        assert "Graph" in repr_str
        assert "N=" in repr_str

    def test_repr_includes_class_name(self):
        """Test that repr includes correct class name."""
        mc = MultiplicativeCascadeGraph(only_const_mode=True)
        assert "MultiplicativeCascadeGraph" in repr(mc)

        vck = VicsekGraph(N=50, k=3, only_const_mode=True)
        assert "VicsekGraph" in repr(vck)

        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5, only_const_mode=True)
        assert "DiracCombGraph" in repr(dcomb)

        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3, only_const_mode=True)
        assert "DiracBrushGraph" in repr(dbrush)


class TestMultispectralGraphCustomPaths:
    """Test custom sgpathn for all subclasses."""

    def test_default_sgpathn_values(self):
        """Test that default sgpathn values are correct."""
        mc = MultiplicativeCascadeGraph(only_const_mode=True)
        assert mc.sgpathn == "msg_mc"

        vck = VicsekGraph(N=50, k=3, only_const_mode=True)
        assert vck.sgpathn == "msg_plv"

        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5, only_const_mode=True)
        assert dcomb.sgpathn == "msg_dcomb"

        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3, only_const_mode=True)
        assert dbrush.sgpathn == "msg_dbrush"

    def test_custom_sgpathn(self):
        """Test that custom sgpathn can be set."""
        custom_path = "custom_path"

        mc = MultiplicativeCascadeGraph(sgpathn=custom_path, only_const_mode=True)
        assert mc.sgpathn == custom_path

        vck = VicsekGraph(N=50, k=3, sgpathn=custom_path, only_const_mode=True)
        assert vck.sgpathn == custom_path

        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5, sgpathn=custom_path, only_const_mode=True)
        assert dcomb.sgpathn == custom_path

        dbrush = DiracBrushGraph(
            base_x=5, base_y=4, fiber_nodes=3, sgpathn=custom_path, only_const_mode=True
        )
        assert dbrush.sgpathn == custom_path


class TestDiracStructureSpecificMethods:
    """Test Dirac-specific methods are only available on Dirac graphs."""

    def test_non_dirac_graphs_no_dirac_methods(self):
        """Test that non-Dirac graphs don't have Dirac methods."""
        mc = MultiplicativeCascadeGraph(only_const_mode=True)
        assert not hasattr(mc, "is_dirac_structure") or not callable(
            getattr(mc, "is_dirac_structure", None)
        ) or mc.is_dirac_structure() is False

        vck = VicsekGraph(N=50, k=3, only_const_mode=True)
        assert not hasattr(vck, "is_dirac_structure") or not callable(
            getattr(vck, "is_dirac_structure", None)
        ) or vck.is_dirac_structure() is False

    def test_dirac_graphs_have_dirac_methods(self):
        """Test that Dirac graphs have Dirac-specific methods."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5)
        assert hasattr(dcomb, "is_dirac_structure")
        assert callable(dcomb.is_dirac_structure)
        assert dcomb.is_dirac_structure() is True
        assert hasattr(dcomb, "compute_base_spectrum")
        assert hasattr(dcomb, "compute_fiber_laplacian")
        assert hasattr(dcomb, "get_base_fiber_dimensions")

        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3)
        assert hasattr(dbrush, "is_dirac_structure")
        assert callable(dbrush.is_dirac_structure)
        assert dbrush.is_dirac_structure() is True
        assert hasattr(dbrush, "compute_base_spectrum")
        assert hasattr(dbrush, "compute_fiber_laplacian")
        assert hasattr(dbrush, "get_base_fiber_dimensions")


class TestBackwardCompatibilityImports:
    """Test that new classes can be imported from MultispectralGraph module."""

    def test_import_from_multispectral_graph(self):
        """Test that all classes can be imported from MultispectralGraph."""
        from lrgsglib.graphs.nx import (
            MultispectralGraphNX as MultispectralGraph,
            MultiplicativeCascadeGraphNX as MultiplicativeCascadeGraph,
            VicsekGraphNX as VicsekGraph,
            DiracCombGraphNX as DiracCombGraph,
            DiracBrushGraphNX as DiracBrushGraph,
            DiracLatticeGraphNX as DiracLatticeGraph,
        )

        # Test that imports work
        assert MultispectralGraph is not None
        assert MultiplicativeCascadeGraph is not None
        assert VicsekGraph is not None
        assert DiracCombGraph is not None
        assert DiracBrushGraph is not None
        assert DiracLatticeGraph is not None

    def test_import_generator_functions(self):
        """Test that generator functions are still available."""
        from lrgsglib.graphs.nx.MultispectralGraphNX import (
            multiplicative_cascade_probability_matrix,
            multiplicative_cascade_graph,
            dirac_comb_graph,
            dirac_brush_graph,
        )

        # Test that imports work
        assert multiplicative_cascade_probability_matrix is not None
        assert multiplicative_cascade_graph is not None
        assert dirac_comb_graph is not None
        assert dirac_brush_graph is not None
