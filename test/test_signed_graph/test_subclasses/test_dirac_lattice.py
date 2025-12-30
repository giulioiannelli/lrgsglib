"""
Unit tests for DiracLatticeGraph, DiracCombGraph, and DiracBrushGraph classes.
"""
import pytest
import networkx as nx
import numpy as np

from lrgsglib.nx_patches.DiracLattice import (
    DiracLatticeGraph,
    DiracCombGraph,
    DiracBrushGraph,
)
from lrgsglib.nx_patches.SignedGraph import SignedGraph


class TestDiracCombGraphCreation:
    """Test DiracCombGraph creation and initialization."""

    def test_default_parameters(self):
        """Test creation with default parameters."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5, only_const_mode=True)
        assert isinstance(dcomb, DiracCombGraph)
        assert isinstance(dcomb, DiracLatticeGraph)
        assert isinstance(dcomb, SignedGraph)
        assert dcomb.only_const_mode is True

    def test_custom_parameters(self):
        """Test creation with custom parameters."""
        dcomb = DiracCombGraph(
            base_nodes=20, fiber_nodes=10, base_type="cycle", periodic=False, only_const_mode=True
        )
        assert dcomb.base_nodes == 20
        assert dcomb.fiber_nodes == 10
        assert dcomb.base_type == "cycle"
        assert dcomb.periodic is False

    def test_base_type_line(self):
        """Test with line base type."""
        dcomb = DiracCombGraph(
            base_nodes=15, fiber_nodes=7, base_type="line", only_const_mode=True
        )
        assert dcomb.base_type == "line"

    def test_base_type_cycle(self):
        """Test with cycle base type."""
        dcomb = DiracCombGraph(
            base_nodes=15, fiber_nodes=7, base_type="cycle", only_const_mode=True
        )
        assert dcomb.base_type == "cycle"

    def test_periodic_true(self):
        """Test periodic boundary conditions."""
        dcomb = DiracCombGraph(
            base_nodes=15, fiber_nodes=7, periodic=True, only_const_mode=True
        )
        assert dcomb.periodic is True

    def test_periodic_false(self):
        """Test non-periodic boundary conditions."""
        dcomb = DiracCombGraph(
            base_nodes=15, fiber_nodes=7, periodic=False, only_const_mode=True
        )
        assert dcomb.periodic is False


class TestDiracCombGraphGeneration:
    """Test DiracCombGraph generation."""

    def test_graph_generation(self):
        """Test that graph is generated correctly."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5)
        assert isinstance(dcomb.G, nx.Graph)
        assert dcomb.N > 0
        assert dcomb.dirac_structure is not None

    def test_graph_not_generated_in_const_mode(self):
        """Test that graph is not generated in only_const_mode."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5, only_const_mode=True)
        assert dcomb.G.number_of_nodes() == 0
        assert dcomb.dirac_structure is None

    def test_dirac_structure_metadata(self):
        """Test that dirac_structure contains correct metadata."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5)
        assert "base_nodes" in dcomb.dirac_structure
        assert "fiber_nodes" in dcomb.dirac_structure
        assert "base_graph" in dcomb.dirac_structure
        assert "fiber_graph" in dcomb.dirac_structure

    def test_base_and_fiber_graphs_set(self):
        """Test that base_graph and fiber_graph are set."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5)
        assert dcomb.base_graph is not None
        assert dcomb.fiber_graph is not None
        assert isinstance(dcomb.base_graph, nx.Graph)
        assert isinstance(dcomb.fiber_graph, nx.Graph)

    def test_H_attribute_set(self):
        """Test that H attribute is set after generation."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5)
        assert dcomb.H is not None
        assert isinstance(dcomb.H, nx.Graph)


class TestDiracCombGraphNodeCount:
    """Test DiracCombGraph node count calculations."""

    def test_expected_num_nodes_const_mode(self):
        """Test that expected_num_nodes returns 0 in const mode."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5, only_const_mode=True)
        assert dcomb.get_expected_num_nodes() == 0

    def test_expected_num_nodes_formula(self):
        """Test expected node count calculation."""
        base_nodes = 10
        fiber_nodes = 5
        dcomb = DiracCombGraph(base_nodes=base_nodes, fiber_nodes=fiber_nodes)
        # Formula: base_nodes * (1 + fiber_nodes)
        expected = base_nodes * (1 + fiber_nodes)
        assert dcomb.get_expected_num_nodes() == expected

    def test_actual_equals_expected_nodes(self):
        """Test that actual node count equals expected."""
        base_nodes = 12
        fiber_nodes = 8
        dcomb = DiracCombGraph(base_nodes=base_nodes, fiber_nodes=fiber_nodes)
        # Formula: base_nodes * (1 + fiber_nodes)
        expected = base_nodes * (1 + fiber_nodes)
        assert dcomb.N == expected


class TestDiracBrushGraphCreation:
    """Test DiracBrushGraph creation and initialization."""

    def test_default_parameters(self):
        """Test creation with default parameters."""
        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3, only_const_mode=True)
        assert isinstance(dbrush, DiracBrushGraph)
        assert isinstance(dbrush, DiracLatticeGraph)
        assert isinstance(dbrush, SignedGraph)
        assert dbrush.only_const_mode is True

    def test_custom_parameters(self):
        """Test creation with custom parameters."""
        dbrush = DiracBrushGraph(
            base_x=8, base_y=6, fiber_nodes=4, periodic=False, only_const_mode=True
        )
        assert dbrush.base_x == 8
        assert dbrush.base_y == 6
        assert dbrush.fiber_nodes == 4
        assert dbrush.periodic is False

    def test_periodic_true(self):
        """Test periodic boundary conditions."""
        dbrush = DiracBrushGraph(
            base_x=5, base_y=4, fiber_nodes=3, periodic=True, only_const_mode=True
        )
        assert dbrush.periodic is True

    def test_periodic_false(self):
        """Test non-periodic boundary conditions."""
        dbrush = DiracBrushGraph(
            base_x=5, base_y=4, fiber_nodes=3, periodic=False, only_const_mode=True
        )
        assert dbrush.periodic is False


class TestDiracBrushGraphGeneration:
    """Test DiracBrushGraph generation."""

    def test_graph_generation(self):
        """Test that graph is generated correctly."""
        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3)
        assert isinstance(dbrush.G, nx.Graph)
        assert dbrush.N > 0
        assert dbrush.dirac_structure is not None

    def test_graph_not_generated_in_const_mode(self):
        """Test that graph is not generated in only_const_mode."""
        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3, only_const_mode=True)
        assert dbrush.G.number_of_nodes() == 0
        assert dbrush.dirac_structure is None

    def test_dirac_structure_metadata(self):
        """Test that dirac_structure contains correct metadata."""
        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3)
        assert "base_nodes" in dbrush.dirac_structure
        assert "fiber_nodes" in dbrush.dirac_structure
        assert "base_graph" in dbrush.dirac_structure
        assert "fiber_graph" in dbrush.dirac_structure

    def test_base_and_fiber_graphs_set(self):
        """Test that base_graph and fiber_graph are set."""
        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3)
        assert dbrush.base_graph is not None
        assert dbrush.fiber_graph is not None
        assert isinstance(dbrush.base_graph, nx.Graph)
        assert isinstance(dbrush.fiber_graph, nx.Graph)


class TestDiracBrushGraphNodeCount:
    """Test DiracBrushGraph node count calculations."""

    def test_expected_num_nodes_const_mode(self):
        """Test that expected_num_nodes returns 0 in const mode."""
        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3, only_const_mode=True)
        assert dbrush.get_expected_num_nodes() == 0

    def test_expected_num_nodes_formula(self):
        """Test expected node count calculation."""
        base_x = 6
        base_y = 5
        fiber_nodes = 4
        dbrush = DiracBrushGraph(base_x=base_x, base_y=base_y, fiber_nodes=fiber_nodes)
        # Formula: base_x * base_y * (1 + fiber_nodes)
        expected = base_x * base_y * (1 + fiber_nodes)
        assert dbrush.get_expected_num_nodes() == expected

    def test_actual_equals_expected_nodes(self):
        """Test that actual node count equals expected."""
        base_x = 7
        base_y = 6
        fiber_nodes = 3
        dbrush = DiracBrushGraph(base_x=base_x, base_y=base_y, fiber_nodes=fiber_nodes)
        # Formula: base_x * base_y * (1 + fiber_nodes)
        expected = base_x * base_y * (1 + fiber_nodes)
        assert dbrush.N == expected


class TestDiracLatticeGraphMethods:
    """Test DiracLatticeGraph specific methods."""

    def test_is_dirac_structure_comb(self):
        """Test that is_dirac_structure returns True for comb."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5)
        assert dcomb.is_dirac_structure() is True

    def test_is_dirac_structure_brush(self):
        """Test that is_dirac_structure returns True for brush."""
        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3)
        assert dbrush.is_dirac_structure() is True

    def test_get_base_fiber_dimensions_comb(self):
        """Test get_base_fiber_dimensions for comb."""
        base_nodes = 12
        fiber_nodes = 6
        dcomb = DiracCombGraph(base_nodes=base_nodes, fiber_nodes=fiber_nodes)
        base_dim, fiber_dim = dcomb.get_base_fiber_dimensions()
        assert base_dim == base_nodes
        assert fiber_dim == fiber_nodes

    def test_get_base_fiber_dimensions_brush(self):
        """Test get_base_fiber_dimensions for brush."""
        base_x = 8
        base_y = 6
        fiber_nodes = 4
        dbrush = DiracBrushGraph(base_x=base_x, base_y=base_y, fiber_nodes=fiber_nodes)
        base_dim, fiber_dim = dbrush.get_base_fiber_dimensions()
        assert base_dim == base_x * base_y
        assert fiber_dim == fiber_nodes


class TestDiracLatticeSpectralMethods:
    """Test spectral computation methods for Dirac structures."""

    def test_compute_base_spectrum_comb(self):
        """Test base spectrum computation for comb."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5)
        base_eigenvalues = dcomb.compute_base_spectrum()
        assert base_eigenvalues is not None
        assert len(base_eigenvalues) == dcomb.base_nodes

    def test_compute_base_spectrum_brush(self):
        """Test base spectrum computation for brush."""
        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3)
        base_eigenvalues = dbrush.compute_base_spectrum()
        assert base_eigenvalues is not None
        assert len(base_eigenvalues) == dbrush.base_x * dbrush.base_y

    def test_compute_fiber_laplacian_comb(self):
        """Test fiber Laplacian computation for comb."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5)
        fiber_lap = dcomb.compute_fiber_laplacian()
        assert fiber_lap is not None
        assert fiber_lap.shape == (dcomb.fiber_nodes, dcomb.fiber_nodes)

    def test_compute_fiber_laplacian_brush(self):
        """Test fiber Laplacian computation for brush."""
        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3)
        fiber_lap = dbrush.compute_fiber_laplacian()
        assert fiber_lap is not None
        assert fiber_lap.shape == (dbrush.fiber_nodes, dbrush.fiber_nodes)

    def test_compute_dirac_spectrum_separated_comb(self):
        """Test separated spectrum computation for comb."""
        dcomb = DiracCombGraph(base_nodes=8, fiber_nodes=4)
        eigenvalues = dcomb.compute_dirac_spectrum_separated()
        assert eigenvalues is not None
        # Separated spectrum has base_nodes * fiber_nodes eigenvalues
        assert len(eigenvalues) == dcomb.base_nodes * dcomb.fiber_nodes

    def test_compute_dirac_spectrum_separated_brush(self):
        """Test separated spectrum computation for brush."""
        dbrush = DiracBrushGraph(base_x=4, base_y=3, fiber_nodes=2)
        eigenvalues = dbrush.compute_dirac_spectrum_separated()
        assert eigenvalues is not None
        # Separated spectrum has (base_x * base_y) * fiber_nodes eigenvalues
        assert len(eigenvalues) == (dbrush.base_x * dbrush.base_y) * dbrush.fiber_nodes

    def test_base_eigenvalues_cached(self):
        """Test that base eigenvalues are cached."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5)
        eigenvalues1 = dcomb.compute_base_spectrum()
        eigenvalues2 = dcomb.compute_base_spectrum()
        assert eigenvalues1 is eigenvalues2  # Same object (cached)

    def test_fiber_laplacian_cached(self):
        """Test that fiber Laplacian is cached."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5)
        lap1 = dcomb.compute_fiber_laplacian()
        lap2 = dcomb.compute_fiber_laplacian()
        assert lap1 is lap2  # Same object (cached)


class TestDiracLatticeIntegration:
    """Test integration with SignedGraph methods."""

    def test_comb_inherits_from_signed_graph(self):
        """Test that DiracCombGraph inherits from SignedGraph."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5, only_const_mode=True)
        assert isinstance(dcomb, SignedGraph)

    def test_brush_inherits_from_signed_graph(self):
        """Test that DiracBrushGraph inherits from SignedGraph."""
        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3, only_const_mode=True)
        assert isinstance(dbrush, SignedGraph)

    def test_laplacian_matrix_computation_comb(self):
        """Test Laplacian matrix computation for comb."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5)
        lap = dcomb.lap
        assert lap is not None
        assert lap.shape[0] == dcomb.N
        assert lap.shape[1] == dcomb.N

    def test_laplacian_matrix_computation_brush(self):
        """Test Laplacian matrix computation for brush."""
        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3)
        lap = dbrush.lap
        assert lap is not None
        assert lap.shape[0] == dbrush.N
        assert lap.shape[1] == dbrush.N

    def test_spectral_computation_comb(self):
        """Test spectral methods work for comb."""
        dcomb = DiracCombGraph(base_nodes=8, fiber_nodes=4)
        dcomb.compute_laplacian_spectrum()
        assert dcomb.eigv is not None
        assert len(dcomb.eigv) == dcomb.N

    def test_spectral_computation_brush(self):
        """Test spectral methods work for brush."""
        dbrush = DiracBrushGraph(base_x=4, base_y=3, fiber_nodes=2)
        dbrush.compute_laplacian_spectrum()
        assert dbrush.eigv is not None
        assert len(dbrush.eigv) == dbrush.N

    def test_sgpathn_attribute_comb(self):
        """Test that sgpathn is set correctly for comb."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5, only_const_mode=True)
        assert dcomb.sgpathn == "msg_dcomb"

    def test_sgpathn_attribute_brush(self):
        """Test that sgpathn is set correctly for brush."""
        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3, only_const_mode=True)
        assert dbrush.sgpathn == "msg_dbrush"


class TestDiracLatticeRepresentation:
    """Test string representations."""

    def test_repr_comb(self):
        """Test __repr__ for DiracCombGraph."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=5)
        repr_str = repr(dcomb)
        assert "DiracCombGraph" in repr_str
        assert "N=" in repr_str

    def test_repr_brush(self):
        """Test __repr__ for DiracBrushGraph."""
        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=3)
        repr_str = repr(dbrush)
        assert "DiracBrushGraph" in repr_str
        assert "N=" in repr_str


class TestDiracLatticeEdgeCases:
    """Test edge cases for Dirac lattice graphs."""

    def test_comb_small_dimensions(self):
        """Test comb with small dimensions."""
        dcomb = DiracCombGraph(base_nodes=3, fiber_nodes=2)
        assert isinstance(dcomb.G, nx.Graph)
        # 3 * (1 + 2) = 9
        assert dcomb.N == 9

    def test_brush_small_dimensions(self):
        """Test brush with small dimensions."""
        dbrush = DiracBrushGraph(base_x=2, base_y=2, fiber_nodes=2)
        assert isinstance(dbrush.G, nx.Graph)
        # 2 * 2 * (1 + 2) = 12
        assert dbrush.N == 12

    def test_comb_single_fiber_node(self):
        """Test comb with single fiber node."""
        dcomb = DiracCombGraph(base_nodes=10, fiber_nodes=1)
        assert isinstance(dcomb.G, nx.Graph)
        # 10 * (1 + 1) = 20
        assert dcomb.N == 20

    def test_brush_single_fiber_node(self):
        """Test brush with single fiber node."""
        dbrush = DiracBrushGraph(base_x=5, base_y=4, fiber_nodes=1)
        assert isinstance(dbrush.G, nx.Graph)
        # 5 * 4 * (1 + 1) = 40
        assert dbrush.N == 40
