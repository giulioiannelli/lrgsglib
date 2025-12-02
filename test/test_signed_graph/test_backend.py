"""Tests for backend abstraction layer (Phase 2).

Tests the pluggable backend system (NumPy/SciPy/CuPy).
"""

import warnings

import numpy as np
import pytest

from lrgsglib.nx_patches.SignedGraph._backend import (
    ArrayBackend,
    Backend,
    BackendManager,
    CupyBackend,
    NumpyBackend,
    ScipyBackend,
)


class TestBackendEnum:
    """Tests for Backend enum."""

    def test_backend_enum_values(self):
        """Test that Backend enum has correct values."""
        assert Backend.NUMPY.value == "numpy"
        assert Backend.SCIPY.value == "scipy"
        assert Backend.CUPY.value == "cupy"

    def test_backend_enum_members(self):
        """Test that all expected backends are present."""
        expected = {"NUMPY", "SCIPY", "CUPY"}
        actual = {b.name for b in Backend}
        assert actual == expected


class TestNumpyBackend:
    """Tests for NumPy backend implementation."""

    @pytest.fixture
    def backend(self):
        """Get NumPy backend instance."""
        return NumpyBackend()

    def test_numpy_backend_name(self, backend):
        """Test backend has correct name."""
        assert backend.name == "numpy"

    def test_numpy_backend_array(self, backend):
        """Test array creation."""
        arr = backend.array([1, 2, 3])
        assert isinstance(arr, np.ndarray)
        np.testing.assert_array_equal(arr, [1, 2, 3])

    def test_numpy_backend_array_with_dtype(self, backend):
        """Test array creation with dtype."""
        arr = backend.array([1.5, 2.5, 3.5], dtype=np.int32)
        assert arr.dtype == np.int32
        np.testing.assert_array_equal(arr, [1, 2, 3])

    def test_numpy_backend_zeros(self, backend):
        """Test zeros creation."""
        arr = backend.zeros((2, 3))
        assert arr.shape == (2, 3)
        np.testing.assert_array_equal(arr, np.zeros((2, 3)))

    def test_numpy_backend_ones(self, backend):
        """Test ones creation."""
        arr = backend.ones((2, 3))
        assert arr.shape == (2, 3)
        np.testing.assert_array_equal(arr, np.ones((2, 3)))

    def test_numpy_backend_eye(self, backend):
        """Test identity matrix creation."""
        arr = backend.eye(3)
        np.testing.assert_array_equal(arr, np.eye(3))

    def test_numpy_backend_dot(self, backend):
        """Test dot product."""
        a = backend.array([1, 2, 3])
        b = backend.array([4, 5, 6])
        result = backend.dot(a, b)
        assert result == 32  # 1*4 + 2*5 + 3*6

    def test_numpy_backend_sum(self, backend):
        """Test sum operation."""
        arr = backend.array([[1, 2], [3, 4]])
        assert backend.sum(arr) == 10
        np.testing.assert_array_equal(backend.sum(arr, axis=0), [4, 6])
        np.testing.assert_array_equal(backend.sum(arr, axis=1), [3, 7])

    def test_numpy_backend_abs(self, backend):
        """Test absolute value."""
        arr = backend.array([-1, -2, 3])
        result = backend.abs(arr)
        np.testing.assert_array_equal(result, [1, 2, 3])

    def test_numpy_backend_sqrt(self, backend):
        """Test square root."""
        arr = backend.array([1, 4, 9])
        result = backend.sqrt(arr)
        np.testing.assert_array_almost_equal(result, [1, 2, 3])

    def test_numpy_backend_eigvalsh(self, backend):
        """Test eigenvalue computation."""
        # Symmetric matrix with known eigenvalues
        mat = backend.array([[2, 1], [1, 2]])
        eigvals = backend.eigvalsh(mat)
        np.testing.assert_array_almost_equal(sorted(eigvals), [1, 3])

    def test_numpy_backend_eigh(self, backend):
        """Test eigenvalue and eigenvector computation."""
        mat = backend.array([[2, 1], [1, 2]])
        eigvals, eigvecs = backend.eigh(mat)

        # Check eigenvalues
        assert len(eigvals) == 2

        # Check that eigenvectors are orthonormal
        assert eigvecs.shape == (2, 2)
        dot_product = np.abs(backend.dot(eigvecs[:, 0], eigvecs[:, 1]))
        assert dot_product < 1e-10  # Nearly orthogonal

    def test_numpy_backend_to_numpy(self, backend):
        """Test conversion to NumPy array."""
        arr = backend.array([1, 2, 3])
        np_arr = backend.to_numpy(arr)
        assert isinstance(np_arr, np.ndarray)
        np.testing.assert_array_equal(np_arr, arr)

    def test_numpy_backend_implements_protocol(self, backend):
        """Test that NumpyBackend implements ArrayBackend protocol."""
        assert isinstance(backend, ArrayBackend)


class TestScipyBackend:
    """Tests for SciPy backend implementation."""

    @pytest.fixture
    def backend(self):
        """Get SciPy backend instance."""
        return ScipyBackend()

    def test_scipy_backend_name(self, backend):
        """Test backend has correct name."""
        assert backend.name == "scipy"

    def test_scipy_backend_eigvalsh(self, backend):
        """Test SciPy eigenvalue computation (uses scipy.linalg)."""
        mat = backend.array([[2, 1], [1, 2]])
        eigvals = backend.eigvalsh(mat)
        np.testing.assert_array_almost_equal(sorted(eigvals), [1, 3])

    def test_scipy_backend_eigh(self, backend):
        """Test SciPy eigenvalue/eigenvector computation."""
        mat = backend.array([[2, 1], [1, 2]])
        eigvals, eigvecs = backend.eigh(mat)

        # Check dimensions
        assert len(eigvals) == 2
        assert eigvecs.shape == (2, 2)

    def test_scipy_backend_implements_protocol(self, backend):
        """Test that ScipyBackend implements ArrayBackend protocol."""
        assert isinstance(backend, ArrayBackend)


class TestCupyBackend:
    """Tests for CuPy backend implementation."""

    def test_cupy_availability_check(self):
        """Test CuPy availability detection."""
        # This test doesn't require CuPy to be installed
        is_available = CupyBackend._check_cupy()
        assert isinstance(is_available, bool)

    @pytest.mark.skipif(
        not CupyBackend._check_cupy(),
        reason="CuPy not installed"
    )
    def test_cupy_backend_name(self):
        """Test CuPy backend name (requires CuPy)."""
        backend = CupyBackend()
        assert backend.name == "cupy"

    @pytest.mark.skipif(
        not CupyBackend._check_cupy(),
        reason="CuPy not installed"
    )
    def test_cupy_backend_array(self):
        """Test CuPy array creation (requires CuPy)."""
        backend = CupyBackend()
        arr = backend.array([1, 2, 3])
        # Convert to numpy for comparison
        np_arr = backend.to_numpy(arr)
        np.testing.assert_array_equal(np_arr, [1, 2, 3])

    def test_cupy_backend_raises_without_cupy(self):
        """Test that CuPy backend raises ImportError when unavailable."""
        if CupyBackend._check_cupy():
            pytest.skip("CuPy is installed, cannot test error path")

        backend = CupyBackend()
        with pytest.raises(ImportError, match="CuPy is not installed"):
            backend.array([1, 2, 3])


class TestBackendManager:
    """Tests for BackendManager."""

    def test_get_numpy_backend(self):
        """Test getting NumPy backend."""
        backend = BackendManager.get_backend(Backend.NUMPY)
        assert backend.name == "numpy"
        assert isinstance(backend, ArrayBackend)

    def test_get_backend_by_string(self):
        """Test getting backend by string name."""
        backend = BackendManager.get_backend("numpy")
        assert backend.name == "numpy"

    def test_get_backend_string_case_insensitive(self):
        """Test that backend string is case-insensitive."""
        backend1 = BackendManager.get_backend("NUMPY")
        backend2 = BackendManager.get_backend("numpy")
        backend3 = BackendManager.get_backend("NumPy")

        # All should return same instance (cached)
        assert backend1 is backend2 is backend3

    def test_get_scipy_backend(self):
        """Test getting SciPy backend."""
        backend = BackendManager.get_backend(Backend.SCIPY)
        assert backend.name == "scipy"

    def test_backend_caching(self):
        """Test that backends are cached."""
        backend1 = BackendManager.get_backend(Backend.NUMPY)
        backend2 = BackendManager.get_backend(Backend.NUMPY)
        assert backend1 is backend2  # Same instance

    def test_invalid_backend_string(self):
        """Test that invalid backend name raises ValueError."""
        with pytest.raises(ValueError, match="Invalid backend"):
            BackendManager.get_backend("invalid_backend")

    def test_cupy_fallback_to_numpy(self):
        """Test automatic fallback to NumPy when CuPy unavailable."""
        if CupyBackend._check_cupy():
            pytest.skip("CuPy is installed, cannot test fallback")

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            backend = BackendManager.get_backend(Backend.CUPY, fallback=True)

            # Should fall back to NumPy
            assert backend.name == "numpy"

            # Should issue warning
            assert len(w) == 1
            assert "CuPy not available" in str(w[0].message)
            assert issubclass(w[0].category, RuntimeWarning)

    def test_cupy_no_fallback_raises(self):
        """Test that CuPy without fallback raises ImportError."""
        if CupyBackend._check_cupy():
            pytest.skip("CuPy is installed, cannot test error")

        with pytest.raises(ImportError, match="CuPy backend requested"):
            BackendManager.get_backend(Backend.CUPY, fallback=False)

    def test_is_available_numpy(self):
        """Test that NumPy backend is always available."""
        assert BackendManager.is_available(Backend.NUMPY)
        assert BackendManager.is_available("numpy")

    def test_is_available_scipy(self):
        """Test that SciPy backend is always available."""
        assert BackendManager.is_available(Backend.SCIPY)
        assert BackendManager.is_available("scipy")

    def test_is_available_cupy(self):
        """Test CuPy availability check."""
        expected = CupyBackend._check_cupy()
        assert BackendManager.is_available(Backend.CUPY) == expected
        assert BackendManager.is_available("cupy") == expected

    def test_is_available_invalid(self):
        """Test that invalid backend name returns False."""
        assert not BackendManager.is_available("invalid")

    def test_list_available(self):
        """Test listing available backends."""
        available = BackendManager.list_available()

        # NumPy and SciPy should always be available
        assert Backend.NUMPY in available
        assert Backend.SCIPY in available

        # CuPy depends on installation
        if CupyBackend._check_cupy():
            assert Backend.CUPY in available
        else:
            assert Backend.CUPY not in available


class TestBackendIntegration:
    """Integration tests for backend system."""

    @pytest.mark.parametrize("backend_name", ["numpy", "scipy"])
    def test_backend_basic_operations(self, backend_name):
        """Test basic operations work across all CPU backends."""
        backend = BackendManager.get_backend(backend_name)

        # Create array
        arr = backend.array([1, 2, 3, 4])

        # Test operations
        assert backend.sum(arr) == 10
        assert backend.sum(backend.abs(backend.array([-1, 2, -3]))) == 6

        # Test matrix operations
        mat = backend.array([[1, 0], [0, 1]])
        np.testing.assert_array_equal(
            backend.to_numpy(mat),
            np.eye(2)
        )

    def test_backend_switching(self):
        """Test switching between backends."""
        # Get NumPy backend
        np_backend = BackendManager.get_backend(Backend.NUMPY)
        arr1 = np_backend.array([1, 2, 3])

        # Get SciPy backend
        sp_backend = BackendManager.get_backend(Backend.SCIPY)
        arr2 = sp_backend.array([1, 2, 3])

        # Both should produce same results
        np.testing.assert_array_equal(
            np_backend.to_numpy(arr1),
            sp_backend.to_numpy(arr2)
        )

    @pytest.mark.skipif(
        not CupyBackend._check_cupy(),
        reason="CuPy not installed"
    )
    def test_cupy_backend_operations(self):
        """Test CuPy backend operations (GPU)."""
        backend = BackendManager.get_backend(Backend.CUPY)

        # Test basic array operations
        arr = backend.array([1, 2, 3, 4])
        assert backend.sum(arr) == 10

        # Test zeros, ones, eye
        zeros = backend.zeros((2, 2))
        assert backend.sum(zeros) == 0

        ones = backend.ones((2, 2))
        assert backend.sum(ones) == 4

        eye = backend.eye(3)
        assert backend.sum(eye) == 3

        # Test mathematical operations
        arr_neg = backend.array([-1, -2, 3])
        abs_arr = backend.abs(arr_neg)
        np.testing.assert_array_equal(
            backend.to_numpy(abs_arr),
            [1, 2, 3]
        )

        sqrt_arr = backend.sqrt(backend.array([1, 4, 9]))
        np.testing.assert_array_almost_equal(
            backend.to_numpy(sqrt_arr),
            [1, 2, 3]
        )

    @pytest.mark.skipif(
        not CupyBackend._check_cupy(),
        reason="CuPy not installed"
    )
    def test_cupy_eigenvalue_computation(self):
        """Test CuPy eigenvalue computation on GPU."""
        backend = BackendManager.get_backend(Backend.CUPY)

        # Symmetric matrix with known eigenvalues
        mat = backend.array([[2.0, 1.0], [1.0, 2.0]])

        # Test eigvalsh
        eigvals = backend.eigvalsh(mat)
        eigvals_np = backend.to_numpy(eigvals)
        np.testing.assert_array_almost_equal(sorted(eigvals_np), [1, 3])

        # Test eigh
        eigvals, eigvecs = backend.eigh(mat)
        eigvals_np = backend.to_numpy(eigvals)
        eigvecs_np = backend.to_numpy(eigvecs)

        assert len(eigvals_np) == 2
        assert eigvecs_np.shape == (2, 2)

    @pytest.mark.skipif(
        not CupyBackend._check_cupy(),
        reason="CuPy not installed"
    )
    def test_numpy_cupy_equivalence(self):
        """Test that NumPy and CuPy backends produce same results."""
        np_backend = BackendManager.get_backend(Backend.NUMPY)
        cp_backend = BackendManager.get_backend(Backend.CUPY)

        # Test array operations
        test_data = [1, 2, 3, 4, 5]

        np_arr = np_backend.array(test_data)
        cp_arr = cp_backend.array(test_data)

        # Compare sums
        assert np_backend.sum(np_arr) == cp_backend.sum(cp_arr)

        # Compare absolute values
        test_data_neg = [-1, 2, -3, 4, -5]
        np_abs = np_backend.abs(np_backend.array(test_data_neg))
        cp_abs = cp_backend.abs(cp_backend.array(test_data_neg))

        np.testing.assert_array_equal(
            np_abs,
            cp_backend.to_numpy(cp_abs)
        )

        # Compare eigenvalues
        mat_data = [[2.0, 1.0], [1.0, 2.0]]
        np_eigvals = np_backend.eigvalsh(np_backend.array(mat_data))
        cp_eigvals = cp_backend.eigvalsh(cp_backend.array(mat_data))

        np.testing.assert_array_almost_equal(
            sorted(np_eigvals),
            sorted(cp_backend.to_numpy(cp_eigvals))
        )

    @pytest.mark.skipif(
        not CupyBackend._check_cupy(),
        reason="CuPy not installed"
    )
    def test_cupy_backend_manager_integration(self):
        """Test BackendManager with CuPy backend."""
        # Get CuPy backend
        backend = BackendManager.get_backend("cupy")
        assert backend.name == "cupy"

        # Verify it's in available backends
        available = BackendManager.list_available()
        assert Backend.CUPY in available

        # Test availability check
        assert BackendManager.is_available(Backend.CUPY)
        assert BackendManager.is_available("cupy")

    @pytest.mark.skipif(
        not CupyBackend._check_cupy(),
        reason="CuPy not installed"
    )
    def test_all_backends_produce_same_results(self):
        """Test that all three backends produce identical results."""
        np_backend = BackendManager.get_backend(Backend.NUMPY)
        sp_backend = BackendManager.get_backend(Backend.SCIPY)
        cp_backend = BackendManager.get_backend(Backend.CUPY)

        # Test data
        mat_data = [[4.0, 2.0], [2.0, 3.0]]

        # Compute eigenvalues with all backends
        np_eigvals = np_backend.eigvalsh(np_backend.array(mat_data))
        sp_eigvals = sp_backend.eigvalsh(sp_backend.array(mat_data))
        cp_eigvals = cp_backend.eigvalsh(cp_backend.array(mat_data))

        # All should match
        np.testing.assert_array_almost_equal(
            sorted(np_eigvals),
            sorted(sp_eigvals),
            decimal=10
        )
        np.testing.assert_array_almost_equal(
            sorted(np_eigvals),
            sorted(cp_backend.to_numpy(cp_eigvals)),
            decimal=10
        )
