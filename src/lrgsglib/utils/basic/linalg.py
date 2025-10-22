import numpy as np
#
from typing import Dict, Iterable, List, Sequence, Tuple, Union
#
from numpy.typing import NDArray
from scipy.spatial.distance import cdist
from scipy.stats import kendalltau, rankdata, spearmanr, wasserstein_distance

__all__ = [
    'basis_random_combination',
    'compute_recon',
    'compute_recon_ultra',
    'compute_mse_from_recon',
    'compute_mse_from_basis',
    'ultrametric_matrix_distance',
    'ultrametric_scaled_distance',
    'ultrametric_rank_correlation',
    'ultrametric_quantile_rmse',
    'compare_ultrametric_trees',
    'ultrametric_distance_permutation_robust',
    'tree_robinson_foulds_distance',
    'tree_cophenetic_correlation',
    'tree_baker_gamma',
    'tree_fowlkes_mallows_index',
    'hierarchical_tree_similarity_suite',
    'ultrametric_multiscale_distance',
    'ultrametric_scale_profile_distance',
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
def ultrametric_matrix_distance(
        D1: NDArray, 
        D2: NDArray, 
        metric: str = 'euclidean'
) -> float:
    """
    Compute the distance between two symmetric matrices using their upper 
    triangular elements.

    This function extracts the upper triangular elements (excluding the 
    diagonal) from both matrices and computes the distance between these 
    flattened vectors using the specified metric. This is particularly useful 
    for comparing distance matrices or correlation matrices where only the 
    upper triangle contains unique information.

    Parameters
    ----------
    D1 : NDArray
        First symmetric matrix of shape (N, N).
    D2 : NDArray
        Second symmetric matrix of shape (N, N). Must have the same shape as D1.
    metric : str, optional
        Distance metric to use for comparing the flattened upper triangular 
        elements. Any metric supported by scipy.spatial.distance.cdist can be 
        used (e.g., 'euclidean', 'manhattan', 'cosine', 'correlation'). 
        Default is 'euclidean'.

    Returns
    -------
    float
        The distance between the two matrices based on their upper triangular 
        elements.

    Notes
    -----
    The function assumes both matrices are symmetric and only uses the upper 
    triangular portion (k=1, excluding diagonal) for comparison. This reduces 
    computational cost and avoids redundant comparisons for symmetric matrices.

    Examples
    --------
    >>> import numpy as np
    >>> D1 = np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
    >>> D2 = np.array([[0, 1.5, 2.1], [1.5, 0, 2.9], [2.1, 2.9, 0]])
    >>> distance = ultrametric_matrix_distance(D1, D2)
    >>> print(f"Euclidean distance: {distance:.3f}")
    """
    # Extract upper triangular indices (excluding diagonal)
    triu_idx = np.triu_indices_from(D1, k=1)
    
    # Flatten upper triangle elements
    v1 = D1[triu_idx]
    v2 = D2[triu_idx]
    
    # Compute distance between flattened vectors
    return cdist([v1], [v2], metric=metric)[0, 0]


def _flatten_ultrametric(M: NDArray) -> NDArray:
    """Return the upper-triangular entries of an ultrametric matrix."""
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise ValueError("Ultrametric matrix must be square.")
    idx = np.triu_indices_from(M, k=1)
    return np.asarray(M, dtype=float)[idx]


def _rescale_ultrametric_vector(
    v: NDArray,
    scale: str = "log",
    eps: float = 1e-12,
    log_base: float = 10.0,
) -> NDArray:
    """Apply scale normalization modes to a flattened ultrametric vector."""
    v = np.asarray(v, dtype=float)

    if scale == "linear":
        return v
    if scale == "log":
        if log_base <= 1:
            raise ValueError("log_base must be greater than 1.")
        return np.log(np.clip(v, eps, None)) / np.log(log_base)
    if scale == "zscore":
        mu = v.mean()
        sigma = v.std()
        return (v - mu) / max(sigma, eps)
    if scale == "rank":
        n = len(v)
        if n == 0:
            return v
        ranks = rankdata(v, method="average") - 1.0
        return ranks / max(n - 1, 1)

    raise ValueError(f"Unknown scale mode '{scale}'.")


def ultrametric_scaled_distance(
    U1: NDArray,
    U2: NDArray,
    metric: str = "euclidean",
    scale: str = "log",
    normalize: bool = True,
    eps: float = 1e-12,
    log_base: float = 10.0,
) -> float:
    """
    Distance between two ultrametric matrices after scale equalisation.

    Let ``U1`` and ``U2`` be ultrametric distance matrices extracted from
    dendrograms on ``n`` leaves. Their unique entries are the upper-triangular
    vectors ``u^(1)``, ``u^(2)`` with components ``u_k = d(i_k, j_k)`` for
    ``k = 1, …, n(n-1)/2``. The ``scale`` argument applies a transformation
    ``T`` to each vector before comparison:

    - ``linear``: ``T(u) = u``.
    - ``log``: ``T(u) = log_{log_base}(max(u, eps))`` to compress orders
      of magnitude while keeping relative ordering.
    - ``zscore``: ``T(u) = (u - μ)/σ`` for location-scale invariance.
    - ``rank``: ``T(u) = (rank(u) - 1)/(m - 1)`` for monotone-invariant scores.

    When ``metric='euclidean'`` and ``normalize=True`` we return the relative
    error

    ``||T(u^(1)) - T(u^(2))||_2 / (||T(u^(1))||_2 + ||T(u^(2))||_2)``,

    which bounds the score in ``[0, 1]`` and prevents the larger tree from
    dominating. For other metrics we defer to
    ``scipy.spatial.distance.cdist`` applied to the transformed vectors.
    """
    if U1.shape != U2.shape:
        raise ValueError("Both ultrametric matrices must have the same shape.")

    v1 = _rescale_ultrametric_vector(_flatten_ultrametric(U1), scale, eps, log_base)
    v2 = _rescale_ultrametric_vector(_flatten_ultrametric(U2), scale, eps, log_base)

    if metric == "euclidean":
        diff = np.linalg.norm(v1 - v2)
        if normalize:
            denom = np.linalg.norm(v1) + np.linalg.norm(v2)
            denom = max(denom, eps)
            diff /= denom
        return float(diff)

    return float(cdist([v1], [v2], metric=metric)[0, 0])


def ultrametric_rank_correlation(
    U1: NDArray,
    U2: NDArray,
    method: str = "spearman",
) -> float:
    """Rank concordance between ultrametric distance profiles.

    The statistic compares the vectors ``u^(1)`` and ``u^(2)`` of unique
    cophenetic distances. A Spearman score evaluates the Pearson correlation of
    the rank-transformed vectors, i.e. the monotone association between the
    two trees. The Kendall option computes Kendall's ``τ`` counting concordant
    minus discordant pairs ``(u_i - u_j)(v_i - v_j)``. Both statistics take
    values in ``[-1, 1]`` and capture topological agreement independently of
    absolute magnitudes.
    """
    if U1.shape != U2.shape:
        raise ValueError("Both ultrametric matrices must have the same shape.")

    v1 = _flatten_ultrametric(U1)
    v2 = _flatten_ultrametric(U2)

    if method == "spearman":
        corr, _ = spearmanr(v1, v2)
    elif method == "kendall":
        corr, _ = kendalltau(v1, v2)
    else:
        raise ValueError("method must be either 'spearman' or 'kendall'.")

    return float(corr)


def ultrametric_quantile_rmse(
    U1: NDArray,
    U2: NDArray,
    quantiles: Sequence[float] = (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0),
    scale: str = "log",
    eps: float = 1e-12,
    log_base: float = 10.0,
) -> float:
    """Root-mean-square error between rescaled ultrametric quantile profiles.

    After transforming each vector ``u`` with ``T`` (see
    :func:`ultrametric_scaled_distance`), we sort the entries and evaluate the
    empirical quantile function ``Q_u(q)`` at the requested ``quantiles``. The
    returned score is

    ``sqrt(mean_q (Q_{u^(1)}(q) - Q_{u^(2)}(q))^2)``,

    which highlights distributional discrepancies: small ``q`` emphasises the
    fine-scale merges, while large ``q`` focuses on coarse merges. Providing a
    dense quantile grid approximates an ``L^2`` distance between the empirical
    cumulative distributions of the transformed ultrametric values.
    """
    if U1.shape != U2.shape:
        raise ValueError("Both ultrametric matrices must have the same shape.")

    if len(quantiles) < 2:
        raise ValueError("At least two quantiles are required.")

    q = np.asarray(quantiles, dtype=float)
    if (q < 0).any() or (q > 1).any():
        raise ValueError("Quantiles must be within [0, 1].")

    q = np.unique(q)
    if q[0] != 0.0 or q[-1] != 1.0:
        q = np.concatenate(([0.0], q, [1.0]))

    v1 = np.sort(_rescale_ultrametric_vector(_flatten_ultrametric(U1), scale, eps, log_base))
    v2 = np.sort(_rescale_ultrametric_vector(_flatten_ultrametric(U2), scale, eps, log_base))

    q1 = np.quantile(v1, q)
    q2 = np.quantile(v2, q)

    return float(np.sqrt(np.mean((q1 - q2) ** 2)))


def compare_ultrametric_trees(
    U1: NDArray,
    U2: NDArray,
    *,
    scale: str = "log",
    value_metric: str = "euclidean",
    normalize: bool = True,
    quantiles: Sequence[float] = (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0),
    eps: float = 1e-12,
    log_base: float = 10.0,
) -> Dict[str, float]:
    """Bundle complementary statistics for comparing two ultrametric trees.

    The dictionary contains:

    - ``scaled_distance`` – the value metric from
      :func:`ultrametric_scaled_distance`, capturing magnitude differences after
      stabilising the scale.
    - ``rank_spearman`` and ``rank_kendall`` – topology-sensitive
      correlations measuring the agreement of merge ordering across all leaf
      pairs.
    - ``quantile_rmse`` – the quantile root-mean-square error from
      :func:`ultrametric_quantile_rmse`, describing the mismatch between the
      multiscale spectra of the trees.
    - ``wasserstein`` – the 1-wasserstein (earth mover's) distance between the
      transformed distance distributions, i.e. ``∫_0^1 |Q_{u^(1)}(q) - Q_{u^(2)}(q)| dq``.

    Combined, these metrics disentangle value discrepancies, ordering agreement,
    and distributional differences so that hierarchies spanning several orders
    of magnitude can be compared without overweighting only the largest merges.
    """
    if U1.shape != U2.shape:
        raise ValueError("Both ultrametric matrices must have the same shape.")

    flat1 = _flatten_ultrametric(U1)
    flat2 = _flatten_ultrametric(U2)

    scaled1 = _rescale_ultrametric_vector(flat1, scale, eps, log_base)
    scaled2 = _rescale_ultrametric_vector(flat2, scale, eps, log_base)

    if value_metric == "euclidean":
        diff = np.linalg.norm(scaled1 - scaled2)
        denom = np.linalg.norm(scaled1) + np.linalg.norm(scaled2)
        if normalize:
            denom = max(denom, eps)
            diff /= denom
        scaled_distance = float(diff)
    else:
        scaled_distance = float(cdist([scaled1], [scaled2], metric=value_metric)[0, 0])

    spearman_corr, _ = spearmanr(flat1, flat2)
    kendall_corr, _ = kendalltau(flat1, flat2)

    quantile_error = ultrametric_quantile_rmse(
        U1,
        U2,
        quantiles=quantiles,
        scale=scale,
        eps=eps,
        log_base=log_base,
    )

    distribution_distance = float(wasserstein_distance(scaled1, scaled2))

    return {
        "scaled_distance": scaled_distance,
        "rank_spearman": float(spearman_corr),
        "rank_kendall": float(kendall_corr),
        "quantile_rmse": quantile_error,
        "wasserstein": distribution_distance,
    }

from scipy.cluster.hierarchy import optimal_leaf_ordering, leaves_list, cophenet
from scipy.spatial.distance import squareform, cdist

def _ultrametric_matrix_from_linkage(Z: NDArray) -> NDArray:
    """
    Cophenetic (ultrametric) matrix from a linkage.
    """
    coph_condensed = cophenet(Z)  # returns condensed distances if Y=None
    U = squareform(coph_condensed)
    np.fill_diagonal(U, 0.0)
    return U

def _order_from_Z_and_D(Z: NDArray, D_condensed: NDArray | None) -> NDArray:
    """
    Leaf order from linkage; if D_condensed is given, use optimal-leaf-ordering (OLO).
    """
    if D_condensed is not None:
        Z = optimal_leaf_ordering(Z, D_condensed)
    return leaves_list(Z)  # indices in [0..N-1]

def _reindex_matrix(M: NDArray, order: NDArray) -> NDArray:
    return M[np.ix_(order, order)]

def ultrametric_distance_permutation_robust(
    Z1: NDArray,
    Z2: NDArray,
    D1_condensed: NDArray | None = None,
    D2_condensed: NDArray | None = None,
    canonical_labels: NDArray | None = None,
    metric: str = "euclidean",
) -> float:
    """
    Compare two hierarchies by the distance between their cophenetic (ultrametric) matrices,
    controlling for leaf permutations via a common leaf order.

    Parameters
    ----------
    Z1, Z2 : linkage matrices (SciPy format, shape (N-1, 4))
        Must refer to the same set of N leaves with indices 0..N-1.
    D1_condensed, D2_condensed : optional condensed distance vectors (length N*(N-1)/2)
        If provided, they are used for optimal-leaf-ordering of Z1 and Z2 respectively.
        If None, the raw dendrogram leaf order is used.
    canonical_labels : optional array of shape (N,)
        If provided, it defines the *common* target order by sorting these labels
        (e.g., actual node names). Otherwise, the common order is taken as the
        OLO order of Z1 (projected onto indices 0..N-1).
    metric : cdist metric (e.g., 'euclidean', 'cosine', 'correlation', ...)

    Returns
    -------
    float : distance between upper triangles of reindexed ultrametric matrices.
    """
    # Build ultrametric (cophenetic) matrices
    U1 = _ultrametric_matrix_from_linkage(Z1)
    U2 = _ultrametric_matrix_from_linkage(Z2)

    N = U1.shape[0]
    assert U2.shape == (N, N), "Both trees must have the same number of leaves."

    # Compute preferred orders (OLO if D provided; else dendrogram order)
    ord1 = _order_from_Z_and_D(Z1, D1_condensed)
    ord2 = _order_from_Z_and_D(Z2, D2_condensed)

    # Define a *common* canonical order of leaves
    if canonical_labels is not None:
        # canonical order = argsort of provided labels
        canon_ord = np.argsort(np.asarray(canonical_labels))
    else:
        # use Z1's order as canonical
        canon_ord = ord1

    # Map ord1/ord2 to the canonical order.
    # Since leaf indices are 0..N-1 consistently, the permutation that takes identity
    # to canon_ord is simply canon_ord itself.
    Pcanon = canon_ord

    # Reindex U1 and U2 to the canonical order.
    # Note: the raw trees may present leaves in ord1/ord2, but U1/U2 are indexed 0..N-1,
    # so to compare fairly we *only* need to place both matrices on the same index order.
    # Using Pcanon ensures identical row/col meaning for both.
    U1r = _reindex_matrix(U1, Pcanon)
    U2r = _reindex_matrix(U2, Pcanon)

    # Extract upper triangles and compute distance
    triu = np.triu_indices(N, k=1)
    v1 = U1r[triu]
    v2 = U2r[triu]
    return float(cdist([v1], [v2], metric=metric)[0, 0])


def tree_robinson_foulds_distance(Z1: NDArray, Z2: NDArray, normalized: bool = True) -> float:
    """
    Compute the Robinson-Foulds distance between two hierarchical trees.
    
    The RF distance counts the number of splits (bipartitions) that differ between
    two trees. This is a purely topological measure that ignores branch lengths.
    
    Parameters
    ----------
    Z1, Z2 : NDArray
        Linkage matrices in SciPy format (shape (N-1, 4))
    normalized : bool, optional
        If True, normalize by the maximum possible RF distance (2*(N-3) for N leaves)
        
    Returns
    -------
    float
        Robinson-Foulds distance, optionally normalized to [0, 1]
    """
    from scipy.cluster.hierarchy import to_tree
    
    def get_splits(Z):
        """Extract all splits (bipartitions) from a linkage matrix."""
        tree = to_tree(Z, rd=False)
        n_leaves = Z.shape[0] + 1
        splits = set()
        
        def traverse(node, leaf_set=None):
            if leaf_set is None:
                leaf_set = set(range(n_leaves))
            
            if node.is_leaf():
                return {node.id}
            
            left_leaves = traverse(node.left, leaf_set)
            right_leaves = traverse(node.right, leaf_set)
            
            # Add this split (bipartition) to our set
            # Store as frozenset of the smaller partition for canonical representation
            if len(left_leaves) <= len(right_leaves):
                splits.add(frozenset(left_leaves))
            else:
                splits.add(frozenset(right_leaves))
                
            return left_leaves | right_leaves
        
        traverse(tree)
        return splits
    
    splits1 = get_splits(Z1)
    splits2 = get_splits(Z2)
    
    # RF distance = number of splits in either tree but not both
    rf_distance = len(splits1.symmetric_difference(splits2))
    
    if normalized:
        n_leaves = Z1.shape[0] + 1
        max_rf = 2 * (n_leaves - 3) if n_leaves > 3 else 0
        rf_distance = rf_distance / max_rf if max_rf > 0 else 0.0
    
    return float(rf_distance)


def tree_cophenetic_correlation(Z1: NDArray, Z2: NDArray) -> float:
    """
    Compute the cophenetic correlation coefficient between two hierarchical trees.
    
    This measures the linear correlation between the cophenetic distances of the
    two trees. High correlation means similar hierarchical structure.
    
    Parameters
    ----------
    Z1, Z2 : NDArray
        Linkage matrices in SciPy format
        
    Returns
    -------
    float
        Cophenetic correlation coefficient [-1, 1]
    """
    from scipy.cluster.hierarchy import cophenet
    from scipy.stats.stats import pearsonr
    
    coph1 = cophenet(Z1)
    coph2 = cophenet(Z2)
    
    corr, _ = pearsonr(coph1, coph2)
    return float(corr)


def tree_baker_gamma(Z1: NDArray, Z2: NDArray) -> float:
    """
    Compute Baker's Gamma correlation between two hierarchical trees.
    
    This measures the rank correlation between the distance matrices used to
    build the trees, providing a measure of how well the relative ordering
    of distances is preserved.
    
    Parameters
    ----------
    Z1, Z2 : NDArray
        Linkage matrices in SciPy format
        
    Returns
    -------
    float
        Baker's Gamma coefficient [-1, 1]
        
    Notes
    -----
    This vectorized implementation is much faster than the naive nested loop
    approach, reducing complexity from O(n⁴) to O(n²) for n leaves.
    """
    from scipy.cluster.hierarchy import cophenet
    
    # Get cophenetic distances
    coph1 = cophenet(Z1)
    coph2 = cophenet(Z2)
    
    # Vectorized computation of all pairwise differences
    # Create difference matrices: diff_matrix[i,j] = coph[i] - coph[j]
    coph1 = np.asarray(coph1)
    coph2 = np.asarray(coph2)
    
    # Broadcasting: coph[:, None] - coph[None, :] creates the full difference matrix
    diff1_matrix = coph1[:, None] - coph1[None, :]
    diff2_matrix = coph2[:, None] - coph2[None, :]
    
    # Compute products of differences
    products = diff1_matrix * diff2_matrix
    
    # Only count upper triangle (i < j) to avoid double counting
    # Use triu with k=1 to exclude diagonal
    upper_tri_mask = np.triu(np.ones_like(products, dtype=bool), k=1)
    products_upper = products[upper_tri_mask]
    
    # Count concordant (positive product) and discordant (negative product)
    concordant = np.sum(products_upper > 0)
    discordant = np.sum(products_upper < 0)
    
    total_pairs = concordant + discordant
    if total_pairs == 0:
        return 0.0
    
    gamma = (concordant - discordant) / total_pairs
    return float(gamma)


def tree_fowlkes_mallows_index(Z1: NDArray, Z2: NDArray, k: int = None) -> float:
    """
    Compute the Fowlkes-Mallows index between two hierarchical trees.
    
    This measures the similarity between flat clusterings obtained by cutting
    the trees at a specific level. If k is not provided, it uses an optimal
    cut based on the tree structure.
    
    Parameters
    ----------
    Z1, Z2 : NDArray
        Linkage matrices in SciPy format
    k : int, optional
        Number of clusters to cut the trees into. If None, uses an optimal cut.
        
    Returns
    -------
    float
        Fowlkes-Mallows index [0, 1]
    """
    from scipy.cluster.hierarchy import fcluster
    
    n_leaves = Z1.shape[0] + 1
    
    if k is None:
        # Use a reasonable default: cut at 1/4 of the maximum height
        k = max(2, n_leaves // 4)
    
    # Get flat clusterings
    clusters1 = fcluster(Z1, k, criterion='maxclust')
    clusters2 = fcluster(Z2, k, criterion='maxclust')
    
    # Count pairs
    n = len(clusters1)
    tp = 0  # True positives: same cluster in both
    fp = 0  # False positives: same in 1, different in 2
    fn = 0  # False negatives: different in 1, same in 2
    
    for i in range(n):
        for j in range(i + 1, n):
            same1 = clusters1[i] == clusters1[j]
            same2 = clusters2[i] == clusters2[j]
            
            if same1 and same2:
                tp += 1
            elif same1 and not same2:
                fp += 1
            elif not same1 and same2:
                fn += 1
    
    if tp + fp == 0 or tp + fn == 0:
        return 0.0
    
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    
    if precision + recall == 0:
        return 0.0
    
    fm_index = np.sqrt(precision * recall)
    return float(fm_index)


def hierarchical_tree_similarity_suite(
    Z1: NDArray, 
    Z2: NDArray,
    D1_condensed: NDArray = None,
    D2_condensed: NDArray = None,
    canonical_labels: NDArray = None
) -> Dict[str, float]:
    """
    Comprehensive suite of similarity measures between two hierarchical trees.
    
    This function computes multiple complementary measures of tree similarity:
    - Topological measures (Robinson-Foulds distance)
    - Distance-based measures (cophenetic correlation, ultrametric distance)
    - Clustering-based measures (Fowlkes-Mallows index)
    - Rank-based measures (Baker's Gamma, Spearman correlation)
    
    Parameters
    ----------
    Z1, Z2 : NDArray
        Linkage matrices in SciPy format
    D1_condensed, D2_condensed : NDArray, optional
        Original distance matrices (condensed form) for permutation-robust comparison
    canonical_labels : NDArray, optional
        Node labels for canonical ordering
        
    Returns
    -------
    Dict[str, float]
        Dictionary containing all similarity measures:
        - 'robinson_foulds': Normalized RF distance [0, 1] (0 = identical topology)
        - 'cophenetic_correlation': Linear correlation of cophenetic distances [-1, 1]
        - 'baker_gamma': Rank correlation measure [-1, 1]  
        - 'fowlkes_mallows': Clustering similarity index [0, 1]
        - 'ultrametric_distance': Permutation-robust ultrametric distance
        - 'ultrametric_scaled_log': Log-scaled ultrametric distance [0, 1]
        - 'spearman_correlation': Spearman correlation of ultrametric distances [-1, 1]
        - 'kendall_correlation': Kendall tau correlation [-1, 1]
    """
    results = {}
    
    # Topological similarity
    results['robinson_foulds'] = tree_robinson_foulds_distance(Z1, Z2, normalized=True)
    
    # Distance-based similarity
    results['cophenetic_correlation'] = tree_cophenetic_correlation(Z1, Z2)
    results['baker_gamma'] = tree_baker_gamma(Z1, Z2)
    
    # Clustering similarity
    results['fowlkes_mallows'] = tree_fowlkes_mallows_index(Z1, Z2)
    
    # Ultrametric matrix comparisons
    if D1_condensed is not None and D2_condensed is not None:
        results['ultrametric_distance'] = ultrametric_distance_permutation_robust(
            Z1, Z2, D1_condensed, D2_condensed, canonical_labels, metric='euclidean'
        )
    else:
        # Fallback without permutation robustness
        U1 = _ultrametric_matrix_from_linkage(Z1)
        U2 = _ultrametric_matrix_from_linkage(Z2)
        results['ultrametric_distance'] = ultrametric_matrix_distance(U1, U2, metric='euclidean')
    
    # Scaled ultrametric comparisons
    U1 = _ultrametric_matrix_from_linkage(Z1)
    U2 = _ultrametric_matrix_from_linkage(Z2)
    
    results['ultrametric_scaled_log'] = ultrametric_scaled_distance(
        U1, U2, metric='euclidean', scale='log', normalize=True
    )
    
    results['spearman_correlation'] = ultrametric_rank_correlation(
        U1, U2, method='spearman'
    )
    
    results['kendall_correlation'] = ultrametric_rank_correlation(
        U1, U2, method='kendall'
    )
    
    return results


def ultrametric_multiscale_distance(
    Z1: NDArray,
    Z2: NDArray,
    D1_condensed: NDArray | None = None,
    D2_condensed: NDArray | None = None,
    canonical_labels: NDArray | None = None,
    scale_weights: str = "log_uniform",
    n_scales: int = 10
) -> Dict[str, float]:
    """
    Compute scale-aware ultrametric distances that properly handle multiple orders of magnitude.
    
    This function addresses the problem where traditional Euclidean distance between ultrametric
    matrices is dominated by large-scale distances, effectively ignoring fine-scale structure.
    It provides multiple distance measures that account for the multi-scale nature of hierarchical
    trees spanning several orders of magnitude.
    
    Parameters
    ----------
    Z1, Z2 : NDArray
        Linkage matrices in SciPy format (shape (N-1, 4))
    D1_condensed, D2_condensed : NDArray, optional
        Original distance matrices (condensed form) for permutation-robust comparison
    canonical_labels : NDArray, optional
        Node labels for canonical ordering
    scale_weights : str, optional
        Weighting scheme for different scales:
        - "uniform": equal weight to all scales
        - "log_uniform": uniform weight in log space (default)
        - "inverse": weight inversely proportional to scale magnitude
    n_scales : int, optional
        Number of scale bins to use for weighted analysis (default: 10)
        
    Returns
    -------
    Dict[str, float]
        Dictionary with multiple distance measures:
        - "raw_euclidean": standard euclidean distance (dominated by large scales)
        - "log_euclidean": euclidean distance in log space
        - "weighted_euclidean": scale-weighted euclidean distance
        - "multiscale_rmse": RMSE across scale bins
    """
    # Get ultrametric matrices
    U1 = _ultrametric_matrix_from_linkage(Z1)
    U2 = _ultrametric_matrix_from_linkage(Z2)
    
    # Handle permutation robustness if needed
    if D1_condensed is not None and D2_condensed is not None:
        order1 = _order_from_Z_and_D(Z1, D1_condensed)
        order2 = _order_from_Z_and_D(Z2, D2_condensed)
        
        if canonical_labels is not None:  
            canon_order = np.argsort(canonical_labels)
            order1 = canon_order
            order2 = canon_order
            
        U1 = _reindex_matrix(U1, order1)
        U2 = _reindex_matrix(U2, order2)
    
    # Extract upper triangular elements
    triu_idx = np.triu_indices_from(U1, k=1)
    v1 = U1[triu_idx]
    v2 = U2[triu_idx]
    
    # Remove zeros and handle log transform
    nonzero_mask = (v1 > 0) & (v2 > 0)
    v1_nz = v1[nonzero_mask]
    v2_nz = v2[nonzero_mask]
    
    results = {}
    
    # 1. Raw Euclidean (current method - dominated by large scales)
    results["raw_euclidean"] = float(np.linalg.norm(v1 - v2))
    
    # 2. Log-space Euclidean (handles orders of magnitude)
    if len(v1_nz) > 0:
        log_v1 = np.log10(v1_nz + 1e-12)
        log_v2 = np.log10(v2_nz + 1e-12)
        results["log_euclidean"] = float(np.linalg.norm(log_v1 - log_v2))
    else:
        results["log_euclidean"] = 0.0
    
    # 3. Scale-weighted Euclidean
    if len(v1_nz) > 0:
        # Create scale bins in log space
        all_vals = np.concatenate([v1_nz, v2_nz])
        log_min, log_max = np.log10(all_vals.min()), np.log10(all_vals.max())
        
        if log_max > log_min:
            scale_bins = np.logspace(log_min, log_max, n_scales + 1)
            
            weighted_dist = 0.0
            total_weight = 0.0
            
            for i in range(n_scales):
                # Find values in this scale bin
                in_bin1 = (v1_nz >= scale_bins[i]) & (v1_nz < scale_bins[i+1])
                in_bin2 = (v2_nz >= scale_bins[i]) & (v2_nz < scale_bins[i+1])
                
                if np.any(in_bin1) and np.any(in_bin2):
                    bin_v1 = v1_nz[in_bin1]
                    bin_v2 = v2_nz[in_bin2]
                    
                    # Make sure we're comparing same number of elements
                    min_len = min(len(bin_v1), len(bin_v2))
                    if min_len > 0:
                        bin_v1_sorted = np.sort(bin_v1)[:min_len]
                        bin_v2_sorted = np.sort(bin_v2)[:min_len]
                        
                        # Scale-dependent weights
                        if scale_weights == "uniform":
                            weight = 1.0
                        elif scale_weights == "log_uniform":
                            weight = 1.0  # uniform in log space
                        elif scale_weights == "inverse":
                            weight = 1.0 / (scale_bins[i] + 1e-12)
                        else:
                            weight = 1.0
                        
                        bin_dist = np.linalg.norm(bin_v1_sorted - bin_v2_sorted)
                        weighted_dist += weight * bin_dist
                        total_weight += weight
            
            results["weighted_euclidean"] = float(weighted_dist / total_weight if total_weight > 0 else 0.0)
        else:
            results["weighted_euclidean"] = results["raw_euclidean"]
    else:
        results["weighted_euclidean"] = 0.0
    
    # 4. Multi-scale RMSE
    if len(v1_nz) > 0 and log_max > log_min:
        scale_rmses = []
        for i in range(n_scales):
            in_bin1 = (v1_nz >= scale_bins[i]) & (v1_nz < scale_bins[i+1])
            in_bin2 = (v2_nz >= scale_bins[i]) & (v2_nz < scale_bins[i+1])
            
            if np.any(in_bin1) and np.any(in_bin2):
                bin_v1 = v1_nz[in_bin1]
                bin_v2 = v2_nz[in_bin2]
                
                min_len = min(len(bin_v1), len(bin_v2))
                if min_len > 0:
                    bin_v1_sorted = np.sort(bin_v1)[:min_len]
                    bin_v2_sorted = np.sort(bin_v2)[:min_len]
                    rmse = np.sqrt(np.mean((bin_v1_sorted - bin_v2_sorted)**2))
                    scale_rmses.append(rmse)
        
        results["multiscale_rmse"] = float(np.mean(scale_rmses) if scale_rmses else 0.0)
    else:
        results["multiscale_rmse"] = 0.0
    
    return results


def ultrametric_scale_profile_distance(
    Z1: NDArray,
    Z2: NDArray,
    D1_condensed: NDArray | None = None,
    D2_condensed: NDArray | None = None,
    canonical_labels: NDArray | None = None,
    n_quantiles: int = 20
) -> Dict[str, float]:
    """
    Compare ultrametric trees by their scale profiles using distribution-based measures.
    
    This function handles multiple orders of magnitude by comparing the distributions
    of ultrametric distances rather than individual distance values. It works in log
    space to give equal importance to all scales.
    
    Parameters
    ----------
    Z1, Z2 : NDArray
        Linkage matrices in SciPy format
    D1_condensed, D2_condensed : NDArray, optional
        Original distance matrices for permutation-robust comparison
    canonical_labels : NDArray, optional
        Node labels for canonical ordering
    n_quantiles : int, optional
        Number of quantiles to use for profile comparison (default: 20)
        
    Returns
    -------
    Dict[str, float]
        Dictionary containing distribution-based similarity measures:
        - 'quantile_rmse': RMSE between quantile profiles in log space
        - 'wasserstein_distance': Earth mover's distance between log distributions
        - 'ks_statistic': Kolmogorov-Smirnov test statistic
    """
    # Get ultrametric matrices (with permutation handling)
    U1 = _ultrametric_matrix_from_linkage(Z1)
    U2 = _ultrametric_matrix_from_linkage(Z2)
    
    if D1_condensed is not None and D2_condensed is not None:
        order1 = _order_from_Z_and_D(Z1, D1_condensed)
        order2 = _order_from_Z_and_D(Z2, D2_condensed)
        
        if canonical_labels is not None:  
            canon_order = np.argsort(canonical_labels)
            order1 = canon_order
            order2 = canon_order
            
        U1 = _reindex_matrix(U1, order1)
        U2 = _reindex_matrix(U2, order2)
    
    # Extract upper triangular elements
    triu_idx = np.triu_indices_from(U1, k=1)
    v1_raw = U1[triu_idx]
    v2_raw = U2[triu_idx]
    
    # Remove zeros for log analysis
    nonzero_mask = (v1_raw > 0) & (v2_raw > 0)
    v1 = v1_raw[nonzero_mask]
    v2 = v2_raw[nonzero_mask]
    
    if len(v1) == 0:
        return {"quantile_rmse": 0.0, "wasserstein_distance": 0.0, "ks_statistic": 0.0}
    
    # Work in log space for multi-scale analysis
    log_v1 = np.log10(v1 + 1e-12)
    log_v2 = np.log10(v2 + 1e-12)
    
    # Quantile-based comparison
    quantiles = np.linspace(0, 1, n_quantiles)
    q1 = np.quantile(log_v1, quantiles)
    q2 = np.quantile(log_v2, quantiles)
    quantile_rmse = float(np.sqrt(np.mean((q1 - q2)**2)))
    
    # Wasserstein distance (optimal transport) in log space
    wasserstein_dist = float(wasserstein_distance(log_v1, log_v2))
    
    # Kolmogorov-Smirnov test statistic
    from scipy.stats import ks_2samp
    ks_stat, _ = ks_2samp(log_v1, log_v2)
    
    return {
        "quantile_rmse": quantile_rmse,
        "wasserstein_distance": wasserstein_dist,
        "ks_statistic": float(ks_stat)
    }
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
