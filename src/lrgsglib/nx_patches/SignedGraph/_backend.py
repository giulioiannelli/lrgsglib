"""Backend abstraction layer for pluggable numerical computing backends.

This module provides a pluggable architecture to support multiple numerical
computing backends (NumPy, SciPy, CuPy) for SignedGraph operations. The backend
system allows users to choose between CPU (NumPy/SciPy) and GPU (CuPy) computation
without changing their code.

Example:
    >>> from lrgsglib.nx_patches.SignedGraph._backend import BackendManager, Backend
    >>> backend = BackendManager.get_backend(Backend.NUMPY)
    >>> arr = backend.array([1, 2, 3])
    >>> print(backend.sum(arr))
    6
"""

from __future__ import annotations

from enum import Enum
from typing import TYPE_CHECKING, Any, Protocol, runtime_checkable

if TYPE_CHECKING:
    from typing import Literal

import numpy as np

__all__ = ["Backend", "ArrayBackend", "BackendManager"]


class Backend(Enum):
    """Enumeration of available numerical computing backends.

    Attributes:
        NUMPY: Pure NumPy backend (CPU, always available)
        SCIPY: SciPy backend with sparse matrix support (CPU)
        CUPY: CuPy backend for GPU acceleration (requires CuPy installation)
    """

    NUMPY = "numpy"
    SCIPY = "scipy"
    CUPY = "cupy"


@runtime_checkable
class ArrayBackend(Protocol):
    """Protocol defining the interface for numerical array backends.

    All backend implementations must provide these methods to ensure
    compatibility with SignedGraph operations.
    """

    name: str

    def array(self, obj: Any, dtype: Any = None) -> Any:
        """Create an array from a Python object.

        Args:
            obj: Input object (list, tuple, ndarray, etc.)
            dtype: Desired data type (optional)

        Returns:
            Array object for this backend
        """
        ...

    def zeros(self, shape: tuple[int, ...], dtype: Any = None) -> Any:
        """Create an array filled with zeros.

        Args:
            shape: Shape of the array
            dtype: Desired data type (optional)

        Returns:
            Zero-filled array
        """
        ...

    def ones(self, shape: tuple[int, ...], dtype: Any = None) -> Any:
        """Create an array filled with ones.

        Args:
            shape: Shape of the array
            dtype: Desired data type (optional)

        Returns:
            One-filled array
        """
        ...

    def eye(self, n: int, dtype: Any = None) -> Any:
        """Create an identity matrix.

        Args:
            n: Size of the matrix (n x n)
            dtype: Desired data type (optional)

        Returns:
            Identity matrix
        """
        ...

    def dot(self, a: Any, b: Any) -> Any:
        """Compute dot product of two arrays.

        Args:
            a: First array
            b: Second array

        Returns:
            Dot product
        """
        ...

    def sum(self, a: Any, axis: int | None = None) -> Any:
        """Sum array elements.

        Args:
            a: Input array
            axis: Axis along which to sum (None = all elements)

        Returns:
            Sum of elements
        """
        ...

    def abs(self, a: Any) -> Any:
        """Compute absolute value element-wise.

        Args:
            a: Input array

        Returns:
            Absolute values
        """
        ...

    def sqrt(self, a: Any) -> Any:
        """Compute square root element-wise.

        Args:
            a: Input array

        Returns:
            Square roots
        """
        ...

    def eigvalsh(self, a: Any) -> Any:
        """Compute eigenvalues of a Hermitian/symmetric matrix.

        Args:
            a: Hermitian/symmetric matrix

        Returns:
            Eigenvalues in ascending order
        """
        ...

    def eigh(self, a: Any) -> tuple[Any, Any]:
        """Compute eigenvalues and eigenvectors of Hermitian/symmetric matrix.

        Args:
            a: Hermitian/symmetric matrix

        Returns:
            Tuple of (eigenvalues, eigenvectors)
        """
        ...

    def eigh_sparse(
        self,
        a: Any,
        k: int = None,
        return_eigenvectors: bool = True,
    ) -> tuple[Any, Any]:
        """Compute eigenvalues/eigenvectors using sparse methods.

        Args:
            a: Sparse or dense Hermitian/symmetric matrix
            k: Number of eigenvalues to compute (None = all possible)
            return_eigenvectors: If False, compute eigenvalues only

        Returns:
            Tuple of (eigenvalues, eigenvectors)
        """
        ...

    def to_numpy(self, a: Any) -> np.ndarray:
        """Convert backend array to NumPy array.

        Args:
            a: Backend array

        Returns:
            NumPy array
        """
        ...


class NumpyBackend:
    """NumPy backend implementation (CPU)."""

    name: str = "numpy"

    @staticmethod
    def array(obj: Any, dtype: Any = None) -> np.ndarray:
        return np.array(obj, dtype=dtype)

    @staticmethod
    def zeros(shape: tuple[int, ...], dtype: Any = None) -> np.ndarray:
        return np.zeros(shape, dtype=dtype)

    @staticmethod
    def ones(shape: tuple[int, ...], dtype: Any = None) -> np.ndarray:
        return np.ones(shape, dtype=dtype)

    @staticmethod
    def eye(n: int, dtype: Any = None) -> np.ndarray:
        return np.eye(n, dtype=dtype)

    @staticmethod
    def dot(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return np.dot(a, b)

    @staticmethod
    def sum(a: np.ndarray, axis: int | None = None) -> np.ndarray | float:
        return np.sum(a, axis=axis)

    @staticmethod
    def abs(a: np.ndarray) -> np.ndarray:
        return np.abs(a)

    @staticmethod
    def sqrt(a: np.ndarray) -> np.ndarray:
        return np.sqrt(a)

    @staticmethod
    def eigvalsh(a: np.ndarray) -> np.ndarray:
        return np.linalg.eigvalsh(a)

    @staticmethod
    def eigh(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        # Use divide-and-conquer driver for full spectrum (faster for N > 1000)
        from scipy import linalg
        return linalg.eigh(a, driver='evd')

    @staticmethod
    def eigh_sparse(
        a: Any,
        k: int = None,
        return_eigenvectors: bool = True,
    ) -> tuple[np.ndarray, np.ndarray | None]:
        """Sparse eigendecomposition using iterative methods."""
        from scipy.sparse.linalg import eigsh
        from scipy.sparse import issparse
        
        # Convert to sparse if needed
        if not issparse(a):
            from scipy.sparse import csr_matrix
            a = csr_matrix(a)
        
        n = a.shape[0]
        if k is None or k >= n - 1:
            # Compute maximum possible (N-2 due to ARPACK limitation)
            k = min(n - 2, n)
        
        # Use shift-invert mode for better convergence
        try:
            result = eigsh(
                a,
                k=k,
                which='SM',
                mode='normal',
                return_eigenvectors=return_eigenvectors,
            )
        except TypeError:
            result = eigsh(a, k=k, which='SM', mode='normal')

        if return_eigenvectors:
            eigvals, eigvecs = result
            return eigvals, eigvecs
        if isinstance(result, tuple):
            result = result[0]
        return result, None

    @staticmethod
    def to_numpy(a: np.ndarray) -> np.ndarray:
        return np.asarray(a)


class ScipyBackend:
    """SciPy backend implementation with sparse matrix support (CPU)."""

    name: str = "scipy"

    @staticmethod
    def array(obj: Any, dtype: Any = None) -> np.ndarray:
        return np.array(obj, dtype=dtype)

    @staticmethod
    def zeros(shape: tuple[int, ...], dtype: Any = None) -> np.ndarray:
        return np.zeros(shape, dtype=dtype)

    @staticmethod
    def ones(shape: tuple[int, ...], dtype: Any = None) -> np.ndarray:
        return np.ones(shape, dtype=dtype)

    @staticmethod
    def eye(n: int, dtype: Any = None) -> np.ndarray:
        return np.eye(n, dtype=dtype)

    @staticmethod
    def dot(a: np.ndarray, b: np.ndarray) -> np.ndarray:
        return np.dot(a, b)

    @staticmethod
    def sum(a: np.ndarray, axis: int | None = None) -> np.ndarray | float:
        return np.sum(a, axis=axis)

    @staticmethod
    def abs(a: np.ndarray) -> np.ndarray:
        return np.abs(a)

    @staticmethod
    def sqrt(a: np.ndarray) -> np.ndarray:
        return np.sqrt(a)

    @staticmethod
    def eigvalsh(a: np.ndarray) -> np.ndarray:
        # Use SciPy for better performance on large matrices
        from scipy import linalg
        return linalg.eigvalsh(a, driver='evd')

    @staticmethod
    def eigh(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        # Use SciPy with divide-and-conquer driver (2x faster for large matrices)
        from scipy import linalg
        return linalg.eigh(a, driver='evd')

    @staticmethod
    def eigh_sparse(
        a: Any,
        k: int = None,
        return_eigenvectors: bool = True,
    ) -> tuple[np.ndarray, np.ndarray | None]:
        """Sparse eigendecomposition using iterative methods."""
        from scipy.sparse.linalg import eigsh
        from scipy.sparse import issparse
        
        # Convert to sparse if needed
        if not issparse(a):
            from scipy.sparse import csr_matrix
            a = csr_matrix(a)
        
        n = a.shape[0]
        if k is None or k >= n - 1:
            # Compute maximum possible (N-2 due to ARPACK limitation)
            k = min(n - 2, n)
        
        # Use shift-invert mode for better convergence
        try:
            result = eigsh(
                a,
                k=k,
                which='SM',
                mode='normal',
                return_eigenvectors=return_eigenvectors,
            )
        except TypeError:
            result = eigsh(a, k=k, which='SM', mode='normal')

        if return_eigenvectors:
            eigvals, eigvecs = result
            return eigvals, eigvecs
        if isinstance(result, tuple):
            result = result[0]
        return result, None

    @staticmethod
    def to_numpy(a: np.ndarray) -> np.ndarray:
        return np.asarray(a)


class CupyBackend:
    """CuPy backend implementation for GPU acceleration.

    Note:
        Requires CuPy to be installed. Falls back to NumPy if unavailable.
    """

    name: str = "cupy"
    _cupy_available: bool | None = None
    _cp: Any = None

    @classmethod
    def _check_cupy(cls) -> bool:
        """Check if CuPy is available and cache result."""
        if cls._cupy_available is None:
            try:
                import cupy as cp
                try:
                    if cp.cuda.runtime.getDeviceCount() < 1:
                        cls._cp = None
                        cls._cupy_available = False
                    else:
                        cls._cp = cp
                        cls._cupy_available = True
                except Exception:
                    cls._cp = None
                    cls._cupy_available = False
            except ImportError:
                cls._cupy_available = False
        return cls._cupy_available

    @classmethod
    def _get_cp(cls) -> Any:
        """Get CuPy module, raising error if unavailable."""
        if not cls._check_cupy():
            raise ImportError(
                "CuPy is not installed or no CUDA device is available. "
                "Install with: pip install cupy-cuda11x (replace cuda11x with "
                "your CUDA version)"
            )
        return cls._cp

    @classmethod
    def array(cls, obj: Any, dtype: Any = None) -> Any:
        cp = cls._get_cp()
        return cp.array(obj, dtype=dtype)

    @classmethod
    def zeros(cls, shape: tuple[int, ...], dtype: Any = None) -> Any:
        cp = cls._get_cp()
        return cp.zeros(shape, dtype=dtype)

    @classmethod
    def ones(cls, shape: tuple[int, ...], dtype: Any = None) -> Any:
        cp = cls._get_cp()
        return cp.ones(shape, dtype=dtype)

    @classmethod
    def eye(cls, n: int, dtype: Any = None) -> Any:
        cp = cls._get_cp()
        return cp.eye(n, dtype=dtype)

    @classmethod
    def dot(cls, a: Any, b: Any) -> Any:
        cp = cls._get_cp()
        return cp.dot(a, b)

    @classmethod
    def sum(cls, a: Any, axis: int | None = None) -> Any:
        cp = cls._get_cp()
        return cp.sum(a, axis=axis)

    @classmethod
    def abs(cls, a: Any) -> Any:
        cp = cls._get_cp()
        return cp.abs(a)

    @classmethod
    def sqrt(cls, a: Any) -> Any:
        cp = cls._get_cp()
        return cp.sqrt(a)

    @classmethod
    def eigvalsh(cls, a: Any) -> Any:
        cp = cls._get_cp()
        # Convert to CuPy array if needed
        a_gpu = cp.asarray(a)
        eigvals = cp.linalg.eigvalsh(a_gpu)
        # Return as NumPy array for consistency
        return cp.asnumpy(eigvals)

    @classmethod
    def eigh(cls, a: Any) -> tuple[Any, Any]:
        cp = cls._get_cp()
        # Convert to CuPy array if needed
        a_gpu = cp.asarray(a)
        eigvals, eigvecs = cp.linalg.eigh(a_gpu)
        # Return as NumPy arrays for consistency
        return cp.asnumpy(eigvals), cp.asnumpy(eigvecs)

    @classmethod
    def eigh_sparse(
        cls,
        a: Any,
        k: int = None,
        return_eigenvectors: bool = True,
    ) -> tuple[Any, Any]:
        """Sparse eigendecomposition using GPU sparse methods when possible."""
        from scipy.sparse import issparse

        n = a.shape[0]
        if k is None or k >= n - 1:
            k = min(n - 2, n)

        # Strategy 1: Try GPU sparse eigsh (CuPy has cupyx.scipy.sparse.linalg.eigsh)
        try:
            cp = cls._get_cp()
            import cupyx.scipy.sparse as cusp
            import cupyx.scipy.sparse.linalg as cusp_linalg

            # Convert to CuPy sparse matrix
            if issparse(a):
                a_gpu = cusp.csr_matrix(a)  # Transfer scipy sparse to GPU
            else:
                # Dense CPU array -> sparse GPU
                from scipy.sparse import csr_matrix
                a_cpu_sparse = csr_matrix(a)
                a_gpu = cusp.csr_matrix(a_cpu_sparse)

            # Use GPU sparse eigsh
            try:
                result = cusp_linalg.eigsh(
                    a_gpu,
                    k=k,
                    which='SA',  # 'SA' = smallest algebraic
                    return_eigenvectors=return_eigenvectors,
                )
            except TypeError:
                result = cusp_linalg.eigsh(a_gpu, k=k, which='SA')

            if return_eigenvectors:
                eigvals_gpu, eigvecs_gpu = result
                eigvals = cp.asnumpy(eigvals_gpu)
                eigvecs = cp.asnumpy(eigvecs_gpu)
                return eigvals, eigvecs

            if isinstance(result, tuple):
                result = result[0]
            eigvals = cp.asnumpy(result)
            return eigvals, None

        except Exception as e:
            import warnings
            warnings.warn(
                f"GPU sparse eigsh failed ({e}), falling back to CPU scipy sparse",
                RuntimeWarning,
                stacklevel=2
            )

            # Strategy 2: Fall back to optimized CPU sparse methods
            from scipy.sparse.linalg import eigsh
            if not issparse(a):
                from scipy.sparse import csr_matrix
                a = csr_matrix(a)

            try:
                result = eigsh(a, k=k, which='SM', return_eigenvectors=return_eigenvectors)
            except TypeError:
                result = eigsh(a, k=k, which='SM')

            if return_eigenvectors:
                eigvals, eigvecs = result
                return eigvals, eigvecs
            if isinstance(result, tuple):
                result = result[0]
            return result, None

    @classmethod
    def to_numpy(cls, a: Any) -> np.ndarray:
        cp = cls._get_cp()
        return cp.asnumpy(a)


class BackendManager:
    """Manager for numerical computing backends.

    Provides a centralized way to get backend instances and handles
    automatic fallback when optional backends are unavailable.

    Example:
        >>> manager = BackendManager()
        >>> numpy_backend = manager.get_backend(Backend.NUMPY)
        >>> arr = numpy_backend.array([1, 2, 3])
    """

    _backends: dict[Backend, type[ArrayBackend]] = {
        Backend.NUMPY: NumpyBackend,
        Backend.SCIPY: ScipyBackend,
        Backend.CUPY: CupyBackend,
    }

    _instances: dict[Backend, ArrayBackend] = {}

    @classmethod
    def get_backend(
        cls,
        backend: Backend | str = Backend.NUMPY,
        fallback: bool = True,
    ) -> ArrayBackend:
        """Get a backend instance by name.

        Args:
            backend: Backend enum or string name ("numpy", "scipy", "cupy")
            fallback: If True, fall back to NumPy if requested backend unavailable

        Returns:
            Backend instance implementing ArrayBackend protocol

        Raises:
            ValueError: If backend name is invalid
            ImportError: If backend unavailable and fallback=False

        Example:
            >>> backend = BackendManager.get_backend(Backend.NUMPY)
            >>> backend = BackendManager.get_backend("scipy")
        """
        # Convert string to enum
        if isinstance(backend, str):
            try:
                backend = Backend(backend.lower())
            except ValueError:
                valid = ", ".join(b.value for b in Backend)
                raise ValueError(
                    f"Invalid backend '{backend}'. Valid options: {valid}"
                )

        # Return cached instance if available
        if backend in cls._instances:
            return cls._instances[backend]

        # Try to create new instance
        backend_class = cls._backends[backend]

        # Special handling for CuPy
        if backend == Backend.CUPY:
            if not CupyBackend._check_cupy():
                if fallback:
                    import warnings
                    warnings.warn(
                        "CuPy not available, falling back to NumPy backend",
                        RuntimeWarning,
                        stacklevel=2,
                    )
                    return cls.get_backend(Backend.NUMPY, fallback=False)
                else:
                    raise ImportError(
                        "CuPy backend requested but CuPy is not installed"
                    )

        # Create and cache instance
        instance = backend_class()
        cls._instances[backend] = instance
        return instance

    @classmethod
    def is_available(cls, backend: Backend | str) -> bool:
        """Check if a backend is available.

        Args:
            backend: Backend enum or string name

        Returns:
            True if backend is available, False otherwise

        Example:
            >>> if BackendManager.is_available(Backend.CUPY):
            ...     backend = BackendManager.get_backend(Backend.CUPY)
        """
        # Convert string to enum
        if isinstance(backend, str):
            try:
                backend = Backend(backend.lower())
            except ValueError:
                return False

        # NumPy and SciPy are always available
        if backend in (Backend.NUMPY, Backend.SCIPY):
            return True

        # Check CuPy
        if backend == Backend.CUPY:
            return CupyBackend._check_cupy()

        return False

    @classmethod
    def list_available(cls) -> list[Backend]:
        """List all available backends.

        Returns:
            List of available Backend enums

        Example:
            >>> available = BackendManager.list_available()
            >>> print(available)
            [<Backend.NUMPY: 'numpy'>, <Backend.SCIPY: 'scipy'>]
        """
        return [b for b in Backend if cls.is_available(b)]
