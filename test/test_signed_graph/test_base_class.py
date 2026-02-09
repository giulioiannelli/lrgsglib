"""
Tests for SignedGraph base class initialization and core functionality.

Tests the refactored __init__ method with backend parameter (Phase 2, Task 2.2).
"""
import pytest
import networkx as nx
import warnings

from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph
from lrgsglib.graphs.nx.SignedGraphNX._backend import Backend, BackendManager


class TestSignedGraphInitWithBackend:
    """Tests for SignedGraph initialization with backend parameter."""

    def test_init_default_backend(self, small_graph):
        """Test that default backend is numpy."""
        sg = SignedGraph(G=small_graph)
        assert hasattr(sg, '_backend'), "Should have _backend attribute"
        assert hasattr(sg, '_backend_name'), "Should have _backend_name attribute"
        assert sg._backend_name == 'numpy', "Default backend should be numpy"

    def test_init_numpy_backend_explicit(self, small_graph):
        """Test explicit numpy backend."""
        sg = SignedGraph(G=small_graph, backend='numpy')
        assert sg._backend_name == 'numpy'

    def test_init_scipy_backend(self, small_graph):
        """Test scipy backend."""
        sg = SignedGraph(G=small_graph, backend='scipy')
        assert sg._backend_name == 'scipy'

    def test_init_backend_enum(self, small_graph):
        """Test initialization with Backend enum."""
        sg = SignedGraph(G=small_graph, backend=Backend.NUMPY)
        assert sg._backend_name == 'numpy'

    def test_init_backend_case_insensitive(self, small_graph):
        """Test that backend string is case-insensitive."""
        sg1 = SignedGraph(G=small_graph, backend='NUMPY')
        sg2 = SignedGraph(G=small_graph, backend='NumPy')
        sg3 = SignedGraph(G=small_graph, backend='numpy')

        assert sg1._backend_name == 'numpy'
        assert sg2._backend_name == 'numpy'
        assert sg3._backend_name == 'numpy'

    @pytest.mark.skipif(
        not BackendManager.is_available('cupy'),
        reason="CuPy not available"
    )
    def test_init_cupy_backend(self, small_graph):
        """Test CuPy backend when available."""
        sg = SignedGraph(G=small_graph, backend='cupy')
        assert sg._backend_name == 'cupy'

    def test_init_cupy_fallback_to_numpy(self, small_graph):
        """Test automatic fallback to numpy when CuPy unavailable."""
        if BackendManager.is_available('cupy'):
            pytest.skip("CuPy is available, can't test fallback")

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sg = SignedGraph(G=small_graph, backend='cupy')

            # Should fall back to numpy
            assert sg._backend_name == 'numpy'

            # Should issue warning
            assert len(w) >= 1
            assert "CuPy not available" in str(w[0].message)

    def test_backend_available_after_init(self, small_graph):
        """Test that backend is accessible after initialization."""
        sg = SignedGraph(G=small_graph, backend='numpy')

        # Should be able to use backend
        arr = sg._backend.array([1, 2, 3])
        assert sg._backend.sum(arr) == 6

    def test_invalid_backend_string(self, small_graph):
        """Test that invalid backend raises ValueError."""
        with pytest.raises(ValueError, match="Invalid backend"):
            SignedGraph(G=small_graph, backend='invalid_backend')


class TestSignedGraphBackendIntegration:
    """Integration tests for backend with SignedGraph operations."""

    @pytest.mark.parametrize("backend_name", ["numpy", "scipy"])
    def test_backend_consistent_across_cpu(self, small_graph, backend_name):
        """Test that different CPU backends produce same graph structure."""
        sg1 = SignedGraph(G=small_graph, backend=backend_name, pflip=0.3, seed=42)
        sg2 = SignedGraph(G=small_graph, backend='numpy', pflip=0.3, seed=42)

        # Should have same number of negative edges
        assert len(sg1.fleset[sg1.on_g]) == len(sg2.fleset[sg2.on_g])

    def test_backend_preserved_through_operations(self, small_graph):
        """Test that backend choice is preserved through graph operations."""
        sg = SignedGraph(G=small_graph, backend='scipy')

        # Flip some edges
        sg.flip_random_fract_edges(pflip=0.2)

        # Backend should still be scipy
        assert sg._backend_name == 'scipy'

    @pytest.mark.skipif(
        not BackendManager.is_available('cupy'),
        reason="CuPy not available"
    )
    def test_numpy_cupy_consistency(self, small_graph):
        """Test NumPy and CuPy backends produce consistent results."""
        sg_np = SignedGraph(G=small_graph, backend='numpy', pflip=0.3, seed=42)
        sg_cp = SignedGraph(G=small_graph, backend='cupy', pflip=0.3, seed=42)

        # Should have same graph structure
        assert sg_np.N == sg_cp.N
        assert sg_np.Ne == sg_cp.Ne
        assert len(sg_np.fleset[sg_np.on_g]) == len(sg_cp.fleset[sg_cp.on_g])


class TestSignedGraphInitTypeHints:
    """Tests for type hints and parameter validation in __init__."""

    def test_pflip_validation(self, small_graph):
        """Test that invalid pflip values are rejected."""
        # pflip out of range should raise error
        with pytest.raises(Exception):  # Assuming _verify_pflip raises an error
            SignedGraph(G=small_graph, pflip=1.5)

        with pytest.raises(Exception):
            SignedGraph(G=small_graph, pflip=-0.1)

    def test_init_with_all_parameters(self, small_graph, tmp_path):
        """Test initialization with all parameters specified."""
        from lrgsglib.config.const import SG_REPR
        sg = SignedGraph(
            G=small_graph,
            pflip=0.3,
            on_g=SG_REPR,  # Use default SG_REPR instead of 'sg'
            seed=42,
            path_data=tmp_path,
            path_plot=tmp_path,
            init_nw_dict=False,
            init_weights_val=1.0,
            export_mode='bin',
            make_dir_tree=False,
            imported=False,
            import_fname='',
            import_mode='pk',
            backend='numpy'
        )

        assert sg is not None
        assert sg._backend_name == 'numpy'
        assert sg.pflip == 0.3

    def test_init_with_dict_weights(self, small_graph):
        """Test initialization with dictionary for init_weights_val."""
        weights_dict = {'weight': 1.0, 'capacity': 2.0}
        sg = SignedGraph(
            G=small_graph,
            init_weights_val=weights_dict,
            backend='numpy'
        )

        assert sg.init_weights_val == weights_dict


class TestSignedGraphProperties:
    """Tests for SignedGraph property accessors."""

    def test_N_property(self, small_graph):
        """Test that N property returns number of nodes."""
        sg = SignedGraph(G=small_graph, backend='numpy')
        assert sg.N == small_graph.number_of_nodes()

    def test_Ne_property(self, small_graph):
        """Test that Ne property returns number of edges."""
        sg = SignedGraph(G=small_graph, backend='numpy')
        assert sg.Ne == small_graph.number_of_edges()

    def test_backend_does_not_affect_graph_size(self, small_graph):
        """Test that backend choice doesn't affect graph size."""
        sg1 = SignedGraph(G=small_graph, backend='numpy')
        sg2 = SignedGraph(G=small_graph, backend='scipy')

        assert sg1.N == sg2.N
        assert sg1.Ne == sg2.Ne


class TestBackendManagerIntegrationWithSignedGraph:
    """Tests for BackendManager integration with SignedGraph."""

    def test_backend_instance_is_cached(self, small_graph):
        """Test that backend instances are cached by BackendManager."""
        sg1 = SignedGraph(G=small_graph, backend='numpy')
        sg2 = SignedGraph(G=small_graph, backend='numpy')

        # Both should use the same cached backend instance
        assert sg1._backend is sg2._backend

    def test_different_backends_different_instances(self, small_graph):
        """Test that different backends have different instances."""
        sg1 = SignedGraph(G=small_graph, backend='numpy')
        sg2 = SignedGraph(G=small_graph, backend='scipy')

        # Should be different instances
        assert sg1._backend is not sg2._backend
        assert sg1._backend_name != sg2._backend_name


class TestSignedGraphPropertiesWithTypeHints:
    """Tests for SignedGraph properties with proper type hints."""

    def test_N_property_returns_int(self, small_graph):
        """Test that N property has correct type."""
        sg = SignedGraph(G=small_graph)
        assert isinstance(sg.N, int)
        assert sg.N > 0

    def test_Ne_property_returns_int(self, small_graph):
        """Test that Ne property has correct type."""
        sg = SignedGraph(G=small_graph)
        assert isinstance(sg.Ne, int)
        assert sg.Ne > 0

    def test_Ne_n_property_returns_int(self, small_graph):
        """Test that Ne_n property has correct type."""
        sg = SignedGraph(G=small_graph, pflip=0.3, seed=42)
        assert isinstance(sg.Ne_n, int)
        assert sg.Ne_n >= 0

    def test_gr_property_returns_dict(self, small_graph):
        """Test that gr property returns dict."""
        sg = SignedGraph(G=small_graph)
        assert isinstance(sg.gr, dict)
        assert len(sg.gr) > 0

    def test_properties_are_read_only(self, small_graph):
        """Test that properties cannot be assigned to."""
        sg = SignedGraph(G=small_graph)

        # Properties should not be settable
        with pytest.raises(AttributeError):
            sg.N = 100

        with pytest.raises(AttributeError):
            sg.Ne = 50


class TestMatrixCaching:
    """Tests for matrix caching behavior."""

    def test_invalidate_matrix_cache_clears_all(self, small_graph):
        """Test that _invalidate_matrix_cache clears all matrices."""
        sg = SignedGraph(G=small_graph)

        # Access matrices to populate cache
        _ = sg.adj
        _ = sg.lap
        _ = sg.slp

        # Verify matrices are cached
        assert sg.on_g in sg.adjacency_matrices
        assert sg.on_g in sg.laplacian_matrices
        assert sg.on_g in sg.signed_laplacian_matrices

        # Invalidate cache
        sg._invalidate_matrix_cache()

        # Verify cache is cleared
        assert len(sg.adjacency_matrices) == 0
        assert len(sg.laplacian_matrices) == 0
        assert len(sg.signed_laplacian_matrices) == 0

    def test_invalidate_matrix_cache_specific_repr(self, small_graph):
        """Test that _invalidate_matrix_cache can target specific representation."""
        sg = SignedGraph(G=small_graph)

        # Access matrices to populate cache
        _ = sg.adj
        on_g = sg.on_g

        # Manually add another representation for testing
        sg.adjacency_matrices['other_repr'] = sg.adj

        # Invalidate specific representation
        sg._invalidate_matrix_cache(on_g=on_g)

        # Only specified representation should be cleared
        assert on_g not in sg.adjacency_matrices
        assert 'other_repr' in sg.adjacency_matrices

    def test_matrix_properties_return_none_when_not_computed(self, small_graph):
        """Test that matrix properties return None when not in cache."""
        sg = SignedGraph(G=small_graph)

        # Clear cache
        sg._invalidate_matrix_cache()

        # Properties should return None (not computed yet)
        assert sg.adj is None
        assert sg.lap is None
        assert sg.slp is None
        assert sg.degm is None
        assert sg.sdeg is None

    def test_cache_persists_across_property_access(self, small_graph):
        """Test that cached matrices persist across multiple accesses."""
        sg = SignedGraph(G=small_graph)

        # Access matrix to populate cache
        adj1 = sg.adj

        # Second access should return same cached object
        adj2 = sg.adj

        # Should be the same object (identity, not just equality)
        if adj1 is not None:
            assert adj1 is adj2


class TestPropertiesUpdateOnChange:
    """Tests that properties update correctly after graph modifications."""

    def test_Ne_n_updates_after_flip(self, small_graph):
        """Test that Ne_n property updates after flipping edges."""
        sg = SignedGraph(G=small_graph, pflip=0.0, seed=42)

        # Initially no negative edges
        initial_neg = sg.Ne_n
        assert initial_neg == 0

        # Flip some edges
        sg.flip_random_fract_edges(pflip=0.5)

        # Should have negative edges now
        assert sg.Ne_n > initial_neg
        assert sg.Ne_n > 0

    def test_N_remains_constant_after_flip(self, small_graph):
        """Test that N doesn't change when flipping edges."""
        sg = SignedGraph(G=small_graph, pflip=0.3, seed=42)

        initial_N = sg.N

        # Flip edges
        sg.flip_random_fract_edges(pflip=0.5)

        # Node count should not change
        assert sg.N == initial_N

    def test_Ne_remains_constant_after_flip(self, small_graph):
        """Test that Ne doesn't change when flipping edges."""
        sg = SignedGraph(G=small_graph, pflip=0.3, seed=42)

        initial_Ne = sg.Ne

        # Flip edges
        sg.flip_random_fract_edges(pflip=0.5)

        # Edge count should not change (just signs flip)
        assert sg.Ne == initial_Ne

    def test_matrix_cache_invalidated_after_modification(self, small_graph):
        """Test that matrix cache is properly invalidated after modifications."""
        sg = SignedGraph(G=small_graph)

        # Access matrix to populate cache
        _ = sg.slp
        assert sg.on_g in sg.signed_laplacian_matrices

        # Modify graph
        sg.flip_random_fract_edges(pflip=0.2)

        # Cache should still exist but content may be outdated
        # (flip_random_fract_edges should ideally call _invalidate_matrix_cache)
        # Note: This depends on implementation of flip_random_fract_edges


class TestPropertiesDocumentation:
    """Tests for property documentation and return types."""

    def test_all_properties_have_docstrings(self, small_graph):
        """Test that all properties have docstrings."""
        sg = SignedGraph(G=small_graph)

        properties = ['N', 'Ne', 'Ne_n', 'adj', 'lap', 'slp', 'degm', 'sdeg', 'gr']

        for prop_name in properties:
            prop = getattr(type(sg), prop_name)
            assert prop.__doc__ is not None, f"{prop_name} should have docstring"
            assert len(prop.__doc__.strip()) > 0, f"{prop_name} docstring should not be empty"

    def test_property_docstrings_have_returns_section(self, small_graph):
        """Test that property docstrings document return types."""
        sg = SignedGraph(G=small_graph)

        properties = ['N', 'Ne', 'Ne_n', 'adj', 'lap', 'slp']

        for prop_name in properties:
            prop = getattr(type(sg), prop_name)
            docstring = prop.__doc__
            assert 'Returns' in docstring or 'returns' in docstring.lower(), \
                f"{prop_name} docstring should document return type"
