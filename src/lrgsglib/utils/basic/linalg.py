import numpy as np
#
from typing import List, Tuple, Union
#
from numpy.typing import NDArray
from scipy.spatial.distance import cdist

__all__ = [
    'basis_random_combination',
    'compute_recon',
    'compute_recon_ultra',
    'compute_mse_from_recon',
    'compute_mse_from_basis',
    'ultrametric_matrix_distance',
    'is_orthonormal',
    'matrix_projection',
    'normalize_array',
    'obtain_coeffs',
    'reconstruct_from_projections',
    'versor',
]
#
def basis_random_combination(
        basis: Union[List[NDArray], NDArray]
) -> Tuple[NDArray, NDArray]:
    """
    Generate a random binary combination of the given basis vectors,
    ensuring that the product involves an odd number (at least 3) of factors.

    Each basis vector is assumed to have entries in {+1, -1}. Instead of
    selecting each coefficient independently, this function chooses a random
    odd number k (with k ≥ 3) and then randomly selects k distinct basis 
    vectors. The resulting vector is computed as the elementwise product of 
    the chosen vectors.

    Parameters
    ----------
    basis : list of numpy.ndarray or numpy.matrix
        A list of basis vectors, each represented as a NumPy array or matrix 
        with entries +1 or -1. The list must contain at least 3 vectors.

    Returns
    -------
    tuple
        A tuple (result, coeffs) where:
          - result is a NumPy array containing the resulting vector.
          - coeffs is a NumPy array of binary coefficients (0 or 1) indicating
            which basis vectors were selected (with an odd count ≥ 3).
    
    Raises
    ------
    ValueError
        If the number of basis vectors is less than 3.
    """
    n = len(basis)
    if n < 3:
        raise ValueError("Basis must contain at least 3 vectors.")

    # Determine possible odd numbers from 3 up to n 
    # (if n is even, max odd is n-1)
    max_odd = n if n % 2 == 1 else n - 1
    possible_counts = np.arange(3, max_odd + 1, 2)
    k = int(np.random.choice(possible_counts))
    
    coeffs = np.zeros(n, dtype=int)
    indices = np.random.choice(n, size=k, replace=False)
    coeffs[indices] = 1

    result = np.ones_like(np.asarray(basis[0]))
    for c, vec in zip(coeffs, basis):
        if c:
            result *= np.asarray(vec)
    return result, coeffs
#
def compute_recon(
        vector: NDArray, 
        basis: List[NDArray], 
        binarize: bool = False
) -> NDArray:
    """
    Compute the reconstruction from a given vector and basis.

    Parameters
    ----------
    vector : NDArray
        The input vector.
    basis : List[NDArray]
        List of basis matrices.

    Returns
    -------
    NDArray
        The sign of the reconstructed matrix.
    """
    projections = np.array(matrix_projection(vector, basis))
    recon = reconstruct_from_projections(projections, basis)
    if binarize:
        recon = np.sign(recon)
    return recon
#
def compute_recon_ultra(
    vec: NDArray,
    B: NDArray,
    mode: str = 'numpy'
) -> NDArray:
    """
    Compute cumulative reconstructions of vec onto prefixes of B.

    Parameters
    ----------
    vec : array-like, shape (N,)
        Input vector.
    B : array-like, shape (K, N)
        Basis matrix with K vectors of length N, pre-normalized.
    mode : {'numpy', 'cupy'}, optional
        Backend to use. 'cupy' requires CuPy installed and GPU.

    Returns
    -------
    recon : ndarray
        Array of shape (K-1, N); row j is reconstruction with B[:j+1].
    """
    # select array module
    try:
        if mode == 'cupy':
            import cupy as cp
            xp = cp
        else:
            xp = np
    except ImportError:
        xp = np

    # convert inputs
    vec_x = xp.asarray(vec)
    B_x = xp.asarray(B)

    # projections and weighted basis
    proj = B_x.dot(vec_x)                # (K,)
    weighted = proj[:, None] * B_x       # (K, N)

    # cumulative reconstructions for prefixes
    recon = xp.cumsum(weighted, axis=0)[:-1]  # (K-1, N)

    return recon
#
def compute_mse_from_recon(recon: NDArray, pattern: NDArray) -> NDArray:
    """
    Compute MSE between a sequence of reconstructions and the original pattern.

    Parameters
    ----------
    recon : ndarray, shape (K-1, N)
        Reconstructions for prefix lengths 1..K-1.
    pattern : array-like, shape (N,)
        Original vector to compare against.

    Returns
    -------
    mse : ndarray, shape (K-1,)
        MSE for each reconstruction step.
    """
    norm_factor = recon[-1].max()
    recon_norm = recon / norm_factor
    diff = recon_norm - pattern[None, :]
    return (diff * diff).mean(axis=1)
#
def compute_mse_from_basis(patterns,
                            basis: NDArray,
                            use_tqdm: bool = False) -> NDArray:
    """
    For each pattern, compute its cumulative reconstructions and MSEs using 
    the provided basis.

    Parameters
    ----------
    patterns : array-like, shape (M, N)
        Collection of input vectors.
    basis : array-like, shape (K, N)
        Basis matrix, pre-normalized rows.
    use_tqdm : bool, optional
        If True, show a progress bar.

    Returns
    -------
    mse_matrix : ndarray, shape (M, K-1)
        mse_matrix[i, j] is the MSE using the first j+1 basis vectors for 
        pattern i.
    """
    iterator = range(len(patterns))
    if use_tqdm:
        from tqdm import tqdm
        iterator = tqdm(iterator, desc='computing MSE')

    M = len(patterns)
    K_minus1 = basis.shape[0] - 1
    mse_matrix = np.empty((M, K_minus1), dtype=np.float32)

    for i in iterator:
        vec = np.asarray(patterns[i], dtype=np.float32)
        recon = compute_recon_ultra(vec, basis)
        mse_matrix[i] = compute_mse_from_recon(recon, vec)

    return mse_matrix
#
def ultrametric_matrix_distance(D1, D2, metric='euclidean'):
    # Flatten upper triangle (excluding diagonal)
    triu_idx = np.triu_indices_find_exactfrom(D1, k=1)
    v1 = D1[triu_idx]
    v2 = D2[triu_idx]
    # Compute distance
    return cdist([v1], [v2], metric=metric)[0, 0]
#
def is_orthonormal(basis: NDArray, axis: int = 0) -> bool:
    """
    Check if a set of basis vectors is orthonormal along a specified axis.

    Parameters
    ----------
    basis : NDArray
        A NumPy array containing the basis vectors. The basis vectors should
        be arranged along the specified axis.
    axis : int, optional
        The axis along which the basis vectors are arranged. Default is 0.

    Returns
    -------
    bool
        True if the basis vectors are orthonormal, False otherwise.

    Notes
    -----
    The function checks orthonormality by computing the dot product of the basis
    vectors and comparing it to the identity matrix. The basis is considered
    orthonormal if:
        basis @ basis.T == I
    where `I` is the identity matrix of appropriate size.
    """
    # Move the specified axis to the first dimension for easier computation
    basis = np.moveaxis(basis, axis, 0)

    # Compute the dot product of the basis vectors
    dot_product = basis @ basis.T

    # Check if the dot product is close to the identity matrix
    identity = np.eye(basis.shape[0])
    return np.allclose(dot_product, identity)
#
def matrix_projection(M: NDArray, basis: List[NDArray]) -> List[float]:
    """
    Compute the projection of a matrix onto a set of basis matrices.

    This function computes the projection of matrix M onto each matrix in the 
    provided basis. The projection onto a basis matrix is defined as the 
    normalized inner product, where the normalization is performed using the 
    Frobenius norm of the basis matrix.

    Parameters
    ----------
    M : NDArray
        2D array representing the matrix to be projected.
    basis : List[NDArray]
        List of 2D arrays representing the basis matrices.

    Returns
    -------
    List[float]
        A list of projection values, one for each basis matrix in the input 
        list.

    Notes
    -----
    The projection for each basis matrix B_i is computed as:
        projection_i = (sum(M * B_i)) / ||B_i||_F,
    where ||B_i||_F denotes the Frobenius norm of B_i.
    """
    projections = []
    for B_i in basis:
        inner_product = np.sum(M * B_i)
        norm_Bi = np.linalg.norm(B_i)
        projection_i = inner_product / norm_Bi
        projections.append(projection_i)
    return projections
#
def normalize_array(array: NDArray, axis: int = None) -> NDArray:
    """
    Normalize a NumPy array along a specified axis.

    Parameters
    ----------
    array : NDArray
        The input array to be normalized.
    axis : int, optional
        The axis along which to normalize the array. If None, the array is 
        flattened and normalized globally. Default is None.

    Returns
    -------
    NDArray
        The normalized array.

    Notes
    -----
    The normalization is performed by dividing the array by its L2 norm along 
    the specified axis. If the axis is None, the entire array is normalized.
    """
    norm = np.linalg.norm(array, axis=axis, keepdims=True)
    return array / norm
#
def obtain_coeffs(basis, vector):
    """
    Obtain the coefficients of the linear combination of the basis vectors that 
    yields the given vector.

    This function assumes that the basis vectors form an invertible set. It 
    converts each basis vector to a NumPy array, constructs the matrix with 
    these vectors as columns, and then solves the linear system
    columns, and then solves the linear system

        B * c = vector

    for the coefficient vector c.

    Parameters
    ----------
    basis : list of numpy.ndarray or numpy.matrix
        A list of basis vectors.
    vector : numpy.ndarray
        The target vector assumed to be a linear combination of the basis 
        vectors.

    Returns
    -------
    numpy.ndarray
        A NumPy array containing the coefficients of the linear combination.

    Raises
    ------
    numpy.linalg.LinAlgError
        If the matrix formed by the basis vectors is singular (i.e., 
        non-invertible).
    """
    basis_arrays = [np.asarray(vec) for vec in basis]
    B = np.column_stack(basis_arrays)
    return np.linalg.solve(B, vector)
#
def reconstruct_from_projections(
        projections: List[float], 
        basis: List[NDArray]
) -> NDArray:
    """
    Reconstruct a matrix from its projections onto a basis.

    This function reconstructs a matrix by summing the product of each 
    projection coefficient with its corresponding basis matrix.

    Parameters
    ----------
    projections : List[float]
        List of projection coefficients for each basis matrix.
    basis : List[NDArray]
        List of basis matrices. All basis matrices must have the same shape.

    Returns
    -------
    NDArray
        The reconstructed matrix obtained by summing the scaled basis matrices.
    """
    reconstructed_matrix = np.zeros_like(basis[0])
    for i, B_i in enumerate(basis):
        reconstructed_matrix += projections[i] * B_i
    return reconstructed_matrix
#
def versor(state_time: NDArray) -> NDArray:
    """
    Compute and return the unit vector (versor) of the input vector.

    Parameters
    ----------
    state_time : NDArray
        The input vector to be normalized.

    Returns
    -------
    NDArray
        The normalized vector (unit vector) in the direction of state_time.

    Raises
    ------
    ValueError
        If the input vector is a zero vector.
    """
    norm = np.linalg.norm(state_time)
    if norm == 0:
        raise ValueError("Zero vector cannot be normalized.")
    return state_time / norm
#