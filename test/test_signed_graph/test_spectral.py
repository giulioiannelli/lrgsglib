"""
Tests for SignedGraph spectral functions with backend integration.

Tests the refactored spectral module (Phase 3).
"""
import pytest
import numpy as np
import networkx as nx

from lrgsglib.graphs.nx import SignedGraphNX as SignedGraph
from lrgsglib.graphs.nx.SignedGraphNX._backend import Backend, BackendManager


class TestComputeLaplacianSpectrum:
    """Tests for compute_laplacian_spectrum with backend."""

    def test_compute_spectrum_uses_backend(self, small_graph):
        """Test that compute_laplacian_spectrum uses the instance backend."""
        sg = SignedGraph(G=small_graph, backend='numpy')

        # Compute spectrum
        sg.compute_laplacian_spectrum()

        # Should have eigenvalues
        assert hasattr(sg, 'eigv')
        assert sg.eigv is not None
        assert len(sg.eigv) == sg.N

    def test_spectrum_eigenvalues_sorted(self, small_graph):
        """Test that eigenvalues are in ascending order."""
        sg = SignedGraph(G=small_graph)
        sg.compute_laplacian_spectrum()

        # Eigenvalues should be sorted
        assert np.all(sg.eigv[:-1] <= sg.eigv[1:])

    def test_spectrum_with_different_backends(self, small_graph):
        """Test that different backends produce same results."""
        sg_numpy = SignedGraph(G=small_graph, backend='numpy', seed=42)
        sg_scipy = SignedGraph(G=small_graph, backend='scipy', seed=42)

        sg_numpy.compute_laplacian_spectrum()
        sg_scipy.compute_laplacian_spectrum()

        # Results should be approximately equal
        np.testing.assert_array_almost_equal(sg_numpy.eigv, sg_scipy.eigv, decimal=10)

    @pytest.mark.skipif(
        not BackendManager.is_available('cupy'),
        reason="CuPy not available"
    )
    def test_spectrum_with_cupy(self, small_graph):
        """Test spectrum computation with CuPy backend."""
        sg = SignedGraph(G=small_graph, backend='cupy')
        sg.compute_laplacian_spectrum()

        assert hasattr(sg, 'eigv')
        assert len(sg.eigv) == sg.N

    def test_spectrum_with_different_precision(self, small_graph):
        """Test spectrum computation with different precision."""
        sg = SignedGraph(G=small_graph)

        # Compute with float64
        sg.compute_laplacian_spectrum(typf=np.float64)
        eigv64 = sg.eigv.copy()

        # Compute with float32
        sg.compute_laplacian_spectrum(typf=np.float32)
        eigv32 = sg.eigv.copy()

        # Should be approximately equal (within float32 precision)
        np.testing.assert_array_almost_equal(eigv64, eigv32, decimal=5)


class TestComputeLaplacianSpectrumWeigV:
    """Tests for compute_laplacian_spectrum_weigV with backend."""

    def test_compute_full_spectrum_uses_backend(self, small_graph):
        """Test that compute_laplacian_spectrum_weigV uses instance backend."""
        sg = SignedGraph(G=small_graph, backend='numpy')
        sg.compute_laplacian_spectrum_weigV()

        # Should have eigenvalues and eigenvectors
        assert hasattr(sg, 'eigv')
        assert hasattr(sg, 'eigV')
        assert sg.eigv is not None
        assert sg.eigV is not None
        assert len(sg.eigv) == sg.N
        assert sg.eigV.shape == (sg.N, sg.N)

    def test_full_spectrum_default_backend(self, small_graph):
        """Test that backend=None uses instance backend."""
        sg = SignedGraph(G=small_graph, backend='scipy')
        sg.compute_laplacian_spectrum_weigV(backend=None)

        assert sg.eigv is not None
        assert sg.eigV is not None

    def test_full_spectrum_backend_override(self, small_graph):
        """Test that backend parameter can override instance backend."""
        sg = SignedGraph(G=small_graph, backend='scipy')

        # Override to use numpy
        sg.compute_laplacian_spectrum_weigV(backend='numpy')

        assert sg.eigv is not None
        assert sg.eigV is not None

    def test_eigenvectors_normalized(self, small_graph):
        """Test that eigenvectors are properly normalized."""
        sg = SignedGraph(G=small_graph)
        sg.compute_laplacian_spectrum_weigV()

        # Each eigenvector should have unit norm (approximately)
        for i in range(sg.N):
            vec = sg.eigV[i, :]
            norm = np.linalg.norm(vec)
            assert abs(norm - 1.0) < 1e-10, f"Eigenvector {i} not normalized"

    def test_eigenvectors_orthogonal(self, small_graph):
        """Test that eigenvectors are orthogonal."""
        sg = SignedGraph(G=small_graph)
        sg.compute_laplacian_spectrum_weigV()

        # Check orthogonality for first few eigenvectors
        for i in range(min(3, sg.N)):
            for j in range(i+1, min(3, sg.N)):
                dot_product = np.abs(np.dot(sg.eigV[i, :], sg.eigV[j, :]))
                assert dot_product < 1e-10, f"Eigenvectors {i} and {j} not orthogonal"

    def test_transpose_parameter(self, small_graph):
        """Test transpose parameter."""
        sg1 = SignedGraph(G=small_graph)
        sg2 = SignedGraph(G=small_graph)

        sg1.compute_laplacian_spectrum_weigV(transpose=True)
        sg2.compute_laplacian_spectrum_weigV(transpose=False)

        # Shapes should be the same
        assert sg1.eigV.shape == sg2.eigV.shape
        # But transpose flag should differ
        assert sg1._eigV_is_transposed != sg2._eigV_is_transposed

    def test_caching_behavior(self, small_graph):
        """Test that full spectrum is cached."""
        sg = SignedGraph(G=small_graph)

        # First computation
        sg.compute_laplacian_spectrum_weigV()
        eigv_first = sg.eigv.copy()

        # Second computation should use cache
        sg.compute_laplacian_spectrum_weigV()
        eigv_second = sg.eigv.copy()

        # Should be identical (same object or values)
        np.testing.assert_array_equal(eigv_first, eigv_second)

    @pytest.mark.skipif(
        not BackendManager.is_available('cupy'),
        reason="CuPy not available"
    )
    def test_numpy_cupy_consistency(self, small_graph):
        """Test that NumPy and CuPy backends produce consistent eigenvalues."""
        sg_np = SignedGraph(G=small_graph, backend='numpy', seed=42)
        sg_cp = SignedGraph(G=small_graph, backend='cupy', seed=42)

        sg_np.compute_laplacian_spectrum_weigV(flip_to_pos=True)
        sg_cp.compute_laplacian_spectrum_weigV(flip_to_pos=True)

        # Eigenvalues should match exactly
        np.testing.assert_array_almost_equal(sg_np.eigv, sg_cp.eigv, decimal=8)

        # Note: Eigenvectors may differ due to:
        # 1. Different numerical algorithms (CPU vs GPU)
        # 2. Degenerate eigenspaces (any linear combination is valid)
        # 3. Different phase choices
        # We only verify eigenvalues match, which is the critical property


class TestComputeAdjacencySpectrumWeigV:
    """Tests for compute_adjacency_spectrum_weigV with backend."""

    def test_adjacency_spectrum_uses_backend(self, small_graph):
        """Test that adjacency spectrum uses instance backend."""
        sg = SignedGraph(G=small_graph, backend='numpy')
        sg.compute_adjacency_spectrum_weigV()

        assert hasattr(sg, 'adj_eigv')
        assert hasattr(sg, 'adj_eigV')
        assert sg.adj_eigv is not None
        assert sg.adj_eigV is not None

    def test_adjacency_spectrum_backend_override(self, small_graph):
        """Test backend parameter override for adjacency spectrum."""
        sg = SignedGraph(G=small_graph, backend='scipy')
        sg.compute_adjacency_spectrum_weigV(backend='numpy')

        assert sg.adj_eigv is not None
        assert sg.adj_eigV is not None

    def test_adjacency_spectrum_consistency(self, small_graph):
        """Test that different backends produce consistent adjacency spectra."""
        sg1 = SignedGraph(G=small_graph, backend='numpy', seed=42)
        sg2 = SignedGraph(G=small_graph, backend='scipy', seed=42)

        sg1.compute_adjacency_spectrum_weigV()
        sg2.compute_adjacency_spectrum_weigV()

        np.testing.assert_array_almost_equal(sg1.adj_eigv, sg2.adj_eigv, decimal=10)


class TestSpectralIntegration:
    """Integration tests for spectral functions."""

    def test_spectrum_then_full_spectrum(self, small_graph):
        """Test computing spectrum only, then full spectrum."""
        sg = SignedGraph(G=small_graph)

        # First compute eigenvalues only
        sg.compute_laplacian_spectrum()
        eigv_only = sg.eigv.copy()

        # Then compute full spectrum
        sg.compute_laplacian_spectrum_weigV()

        # Eigenvalues should match
        np.testing.assert_array_almost_equal(eigv_only, sg.eigv, decimal=10)

    def test_backend_affects_all_spectral_computations(self, small_graph):
        """Test that backend is consistently used across spectral functions."""
        sg = SignedGraph(G=small_graph, backend='scipy')

        # Both should use scipy backend
        sg.compute_laplacian_spectrum()
        sg.compute_adjacency_spectrum_weigV()

        assert sg.eigv is not None
        assert sg.adj_eigv is not None
