import numpy as np
#
from typing import Tuple, Callable, Optional, Dict, Any
from networkx import Graph
#
from numpy.typing import NDArray
from scipy.cluster.hierarchy import cophenet
from scipy.linalg import expm, eigh_tridiagonal
from scipy.sparse import csr_matrix, csr_array
from scipy.sparse.linalg import expm as sparse_expm
from scipy.spatial.distance import squareform
#
from ..basic import dtype_numerical_precision
from .spectral import *
#
__all__ = [
    "extract_ultrametric_matrix",
    "lapl_dists",
    "entropy",
    "compute_entropy_observables_from_eigenvalues",
    "compute_entropy_observables_slq",
    "compute_renyi_observables_from_eigenvalues",
]
#
def _normalize_entropy_profile(entropy_profile: NDArray, entropy_norm: str) -> NDArray:
    norm = entropy_norm.lower()
    if norm == "standard":
        return entropy_profile
    if norm == "complement":
        return 1 - entropy_profile
    raise ValueError("entropy_norm must be 'standard' or 'complement'.")


def extract_ultrametric_matrix(linkage_matrix, n_nodes):
    """
    Extract the ultrametric distance matrix from a linkage matrix.
    
    Parameters:
    -----------
    linkage_matrix : np.ndarray
        The linkage matrix from hierarchical clustering
    n_nodes : int
        Number of original nodes/observations
        
    Returns:
    --------
    np.ndarray
        The ultrametric distance matrix (symmetric, n_nodes x n_nodes)
    """
    # Compute cophenetic distances (these are the ultrametric distances)
    cophenetic_dists = cophenet(linkage_matrix)
    
    # Convert from condensed form to square matrix
    ultrametric_matrix = squareform(cophenetic_dists)
    
    return ultrametric_matrix
#
def lapl_dists(L, tau: float = 1e-2, is_signed: bool = False) -> NDArray:
    """
    Compute the pairwise distances between nodes based on the Laplacian spectrum.

    Parameters
    ----------
    L : NDArray or csr_matrix
        The Laplacian matrix of the graph. Can be a dense numpy array or a sparse scipy CSR matrix.
    tau : float, optional
        A scaling parameter for the exponential function (default is 1e-2).
    is_signed : bool, optional
        Whether the graph is signed (default is False).

    Returns
    -------
    NDArray
        A condensed distance matrix representing pairwise distances between nodes.
    """
    if isinstance(L, (csr_matrix, csr_array)):
        # Sparse matrix computation
        num = sparse_expm(-tau * L)
        trace_num = num.diagonal().sum()
        rho = num / trace_num
        Trho = rho.copy()
        Trho.data = 1.0 / rho.data
        Trho = Trho.maximum(Trho.T)
        Trho.setdiag(0)
        if is_signed:
            old_d = squareform(Trho.toarray())
            dists = np.sqrt(1 - old_d/np.max(old_d))
        else:
            dists = squareform(Trho.toarray())
    else:
        # Dense matrix computation
        num = expm((-tau * L))
        rho = num / np.trace(num)
        Trho = np.copy(1.0 / rho)
        Trho = np.maximum(Trho, Trho.T)
        np.fill_diagonal(Trho, 0)
        if is_signed:
            old_d = squareform(Trho)
            dists = np.sqrt(np.max(old_d) - old_d)
        else:
            dists = squareform(Trho)
    return dists


def entropy(
    G: Optional[Graph] = None,
    steps: int = 600,
    t1: int = -2,
    t2: int = 5,
    wTresh: float = dtype_numerical_precision(np.float64),
    func_lspectrum: Optional[Callable] = get_graph_lspectrum,
    eigenvalues: Optional[NDArray] = None,
    num_nodes: Optional[int] = None,
    entropy_norm: str = "complement",
    specific_heat_scale: str = "logN",
    typf: type = np.float64,
) -> Tuple[NDArray, NDArray, NDArray, NDArray]:
    """
    Compute the entropy, its derivative, variance, and time steps for a graph.

    Parameters
    ----------
    G : networkx.Graph, optional
        The input graph. Required if ``eigenvalues`` is not provided.
    steps : int, optional
        Number of time steps for the computation (default is 600).
    t1 : int, optional
        Logarithmic start time (default is -2).
    t2 : int, optional
        Logarithmic end time (default is 5).
    wTresh : float, optional
        Threshold for filtering eigenvalues (default is machine precision for
        float64). Use None to infer from ``typf``.
    func_lspectrum : Callable, optional
        Function to compute the Laplacian spectrum of the graph (default is
        `get_graph_lspectrum`). Required if ``eigenvalues`` is not provided.
    eigenvalues : NDArray, optional
        Precomputed Laplacian spectrum. If provided, ``G`` is optional.
    num_nodes : int, optional
        Number of nodes for entropy normalization. Defaults to ``len(eigenvalues)``
        when ``eigenvalues`` is provided.
    entropy_norm : str, optional
        Entropy normalization mode: "standard" or "complement".
    specific_heat_scale : str, optional
        Scaling for the entropy derivative: "logN" or "none".
    typf : type, optional
        Floating-point type used for computations (default np.float64).

    Returns
    -------
    Tuple[NDArray, NDArray, NDArray, NDArray]
        - Normalized entropy values.
        - Derivative of entropy with respect to time.
        - Variance of the Laplacian spectrum.
        - Time steps used for the computation.
    """
    if wTresh is None:
        wTresh = dtype_numerical_precision(typf)

    if eigenvalues is None:
        if G is None:
            raise ValueError("Provide G or eigenvalues to compute entropy.")
        if func_lspectrum is None:
            raise ValueError("func_lspectrum must be provided when G is used.")
        N = G.number_of_nodes()
        _, w = func_lspectrum(G)
    else:
        w = eigenvalues
        N = num_nodes if num_nodes is not None else int(np.size(w))

    normalized_entropy, entropy_derivative, variance_profile, time_grid = (
        compute_entropy_observables_from_eigenvalues(
            eigenvalues=w,
            num_nodes=N,
            steps=steps,
            t1=t1,
            t2=t2,
            typf=typf,
            threshold=wTresh,
            entropy_norm=entropy_norm,
            specific_heat_scale=specific_heat_scale,
        )
    )

    return normalized_entropy, entropy_derivative, variance_profile, time_grid


def _lanczos_tridiagonal(
    matvec: Callable[[NDArray], NDArray],
    n: int,
    m: int,
    v0: NDArray,
    typf: type,
) -> Tuple[NDArray, NDArray]:
    v_prev = np.zeros(n, dtype=typf)
    v = v0 / np.linalg.norm(v0)
    alphas = np.zeros(m, dtype=typf)
    betas = np.zeros(m - 1, dtype=typf)
    beta = typf(0)

    for j in range(m):
        w = matvec(v)
        if j > 0:
            w = w - beta * v_prev
        alpha = np.dot(v, w)
        w = w - alpha * v
        beta = np.linalg.norm(w)

        alphas[j] = alpha
        if j < m - 1:
            betas[j] = beta

        if beta == 0:
            return alphas[: j + 1], betas[: j]

        v_prev, v = v, w / beta

    return alphas, betas


def _slq_trace_estimates(
    laplacian: NDArray,
    time_grid: NDArray,
    num_samples: int,
    lanczos_steps: int,
    seed: Optional[int],
    typf: type,
) -> Tuple[NDArray, NDArray, NDArray]:
    n = laplacian.shape[0]
    rng = np.random.default_rng(seed)
    trace_exp = np.zeros_like(time_grid, dtype=typf)
    trace_l_exp = np.zeros_like(time_grid, dtype=typf)
    trace_l2_exp = np.zeros_like(time_grid, dtype=typf)

    matvec = laplacian.dot if hasattr(laplacian, "dot") else lambda x: laplacian @ x

    for _ in range(num_samples):
        r = rng.choice([-1.0, 1.0], size=n).astype(typf)
        norm_sq = np.dot(r, r)
        alphas, betas = _lanczos_tridiagonal(
            matvec=matvec,
            n=n,
            m=lanczos_steps,
            v0=r,
            typf=typf,
        )
        evals, evecs = eigh_tridiagonal(alphas, betas)
        weights = evecs[0, :] ** 2

        exp_terms = np.exp(-np.outer(time_grid, evals))
        w0 = exp_terms @ weights
        w1 = exp_terms @ (weights * evals)
        w2 = exp_terms @ (weights * (evals ** 2))

        trace_exp += norm_sq * w0
        trace_l_exp += norm_sq * w1
        trace_l2_exp += norm_sq * w2

    inv_samples = typf(1.0 / num_samples)
    return trace_exp * inv_samples, trace_l_exp * inv_samples, trace_l2_exp * inv_samples


def compute_entropy_observables_slq(
    laplacian: NDArray,
    num_nodes: int,
    steps: int = 600,
    t1: int = -2,
    t2: int = 5,
    typf: type = np.float64,
    entropy_norm: str = "complement",
    specific_heat_scale: str = "logN",
    lanczos_steps: int = 50,
    num_samples: int = 30,
    seed: Optional[int] = None,
) -> Tuple[NDArray, NDArray, NDArray, NDArray]:
    """
    Approximate entropy observables with stochastic Lanczos quadrature.

    This avoids explicit eigenvalue computation by estimating traces of
    exp(-tau L), L exp(-tau L), and L^2 exp(-tau L) with Hutchinson sampling.
    """
    if laplacian.shape[0] != laplacian.shape[1]:
        raise ValueError("laplacian must be a square matrix.")
    if steps < 1:
        raise ValueError("steps must be at least 1 to build the entropy profile.")
    if num_samples < 1:
        raise ValueError("num_samples must be at least 1.")
    if lanczos_steps < 2:
        raise ValueError("lanczos_steps must be at least 2.")

    laplacian = (
        laplacian.astype(typf)
        if hasattr(laplacian, "astype")
        else np.asarray(laplacian, dtype=typf)
    )

    time_grid = np.logspace(t1, t2, steps, dtype=typf)
    trace_exp, trace_l_exp, trace_l2_exp = _slq_trace_estimates(
        laplacian=laplacian,
        time_grid=time_grid,
        num_samples=num_samples,
        lanczos_steps=lanczos_steps,
        seed=seed,
        typf=typf,
    )

    log_N = np.log(typf(num_nodes)) if num_nodes > 0 else typf(1)
    entropy_profile = np.zeros_like(time_grid, dtype=typf)
    variance_profile = np.zeros_like(time_grid, dtype=typf)

    valid = trace_exp > 0
    avg = np.zeros_like(time_grid, dtype=typf)
    avg2 = np.zeros_like(time_grid, dtype=typf)
    avg[valid] = trace_l_exp[valid] / trace_exp[valid]
    avg2[valid] = trace_l2_exp[valid] / trace_exp[valid]
    variance_profile[valid] = avg2[valid] - avg[valid] ** 2

    with np.errstate(divide="ignore", invalid="ignore"):
        entropy_profile[valid] = (
            time_grid[valid] * avg[valid] + np.log(trace_exp[valid])
        ) / log_N

    normalized_entropy = _normalize_entropy_profile(entropy_profile, entropy_norm)
    entropy_derivative = np.gradient(normalized_entropy, np.log(time_grid))
    scale = specific_heat_scale.lower()
    if scale == "logn":
        entropy_derivative = log_N * entropy_derivative
    elif scale != "none":
        raise ValueError("specific_heat_scale must be 'logN' or 'none'.")

    return normalized_entropy, entropy_derivative, variance_profile, time_grid


def compute_entropy_observables_from_eigenvalues(
    eigenvalues: NDArray,
    num_nodes: int,
    steps: int = 600,
    t1: int = -2,
    t2: int = 5,
    typf: type = np.float64,
    threshold: Optional[float] = None,
    entropy_norm: str = "complement",
    specific_heat_scale: str = "logN",
) -> Tuple[NDArray, NDArray, NDArray, NDArray]:
    """
    Compute entropy-related observables from a Laplacian spectrum.

    Parameters
    ----------
    eigenvalues : NDArray
        Spectrum of the Laplacian (or signed Laplacian) matrix.
    num_nodes : int
        Number of nodes in the graph used to normalize the entropy.
    steps : int, optional
        Number of logarithmic time samples (default 600).
    t1 : int, optional
        Lower exponent for the log-spaced time grid (default -2).
    t2 : int, optional
        Upper exponent for the log-spaced time grid (default 5).
    typf : type, optional
        Floating-point type used for computations (default np.float64).
    threshold : float, optional
        Magnitude threshold below which eigenvalues are treated as zero.
        Defaults to the machine precision of `typf`.
    entropy_norm : str, optional
        Entropy normalization mode: "standard" or "complement".
    specific_heat_scale : str, optional
        Scaling for the entropy derivative: "logN" or "none".

    Returns
    -------
    Tuple[NDArray, NDArray, NDArray, NDArray]
        Normalized entropy profile, entropy derivative, variance profile,
        and the sampled time grid.
    """
    if steps < 1:
        raise ValueError("steps must be at least 1 to build the entropy profile.")

    eigvals = np.asarray(eigenvalues, dtype=typf)
    eps = threshold if threshold is not None else dtype_numerical_precision(typf)
    eigvals = np.where(np.abs(eigvals) > eps, eigvals, typf(0))

    time_grid = np.logspace(t1, t2, steps, dtype=typf)
    entropy_profile = np.zeros_like(time_grid, dtype=typf)
    variance_profile = np.zeros_like(time_grid, dtype=typf)

    log_N = np.log(typf(num_nodes)) if num_nodes > 0 else typf(1)

    for idx, tau in enumerate(time_grid):
        rhoTr = np.exp(-tau * eigvals)
        trace_rho = np.nansum(rhoTr, dtype=typf)
        if trace_rho == 0:
            rho = np.zeros_like(rhoTr)
        else:
            rho = rhoTr / trace_rho

        with np.errstate(divide='ignore', invalid='ignore'):
            entropy_profile[idx] = -np.nansum(rho * np.log(rho), dtype=typf) / log_N

        if trace_rho:
            avgrho = np.nansum(eigvals * rhoTr, dtype=typf) / trace_rho
            av2rho = np.nansum((eigvals ** 2) * rhoTr, dtype=typf) / trace_rho
        else:
            avgrho = typf(0)
            av2rho = typf(0)
        variance_profile[idx] = av2rho - avgrho ** 2

    normalized_entropy = _normalize_entropy_profile(entropy_profile, entropy_norm)

    entropy_derivative = np.gradient(normalized_entropy, np.log(time_grid))
    scale = specific_heat_scale.lower()
    if scale == "logn":
        entropy_derivative = log_N * entropy_derivative
    elif scale != "none":
        raise ValueError("specific_heat_scale must be 'logN' or 'none'.")
    return normalized_entropy, entropy_derivative, variance_profile, time_grid


def compute_renyi_observables_from_eigenvalues(
    eigenvalues: NDArray,
    num_nodes: int,
    q: float,
    steps: int = 600,
    t1: int = -2,
    t2: int = 5,
    typf: type = np.float64,
    threshold: Optional[float] = None,
    tail_fraction: float = 0.1,
    drop_invalid: bool = True,
    entropy_norm: Optional[str] = None,
    legacy_normalization: bool = False,
    specific_heat_scale: str = "none",
) -> Dict[str, Any]:
    """
    Compute Renyi/Tsallis entropy profiles and generalized specific heat.

    Parameters
    ----------
    eigenvalues : NDArray
        Spectrum of the Laplacian (or signed Laplacian) matrix.
    num_nodes : int
        Number of nodes in the graph used to normalize the entropy.
    q : float
        Order of the Renyi/Tsallis entropy (must be > 0).
    steps : int, optional
        Number of logarithmic time samples (default 600).
    t1 : int, optional
        Lower exponent for the log-spaced time grid (default -2).
    t2 : int, optional
        Upper exponent for the log-spaced time grid (default 5).
    typf : type, optional
        Floating-point type used for computations (default np.float64).
    threshold : float, optional
        Magnitude threshold below which eigenvalues are treated as zero.
        Defaults to the machine precision of `typf`.
    tail_fraction : float, optional
        Fraction of the tail used for ds estimate (default 0.1).
    drop_invalid : bool, optional
        If True, discard invalid weights for q != 1 (default True).
    entropy_norm : str, optional
        Entropy normalization mode: "standard" or "complement". If None,
        defaults to "standard" in legacy mode and "complement" otherwise.
    legacy_normalization : bool, optional
        Compatibility mode for older results. When True and ``entropy_norm`` is
        not provided, the entropy profile is returned in "standard" form while
        the specific heat is computed from the complement.
    specific_heat_scale : str, optional
        Scaling for the entropy derivative: "logN" or "none".

    Returns
    -------
    dict
        Dictionary with keys:
        - "q": q value
        - "tau": time grid
        - "renyi_entropy": normalized Renyi entropy profile
        - "tsallis_entropy": normalized Tsallis entropy profile
        - "specific_heat": generalized specific heat profile
        - "ds_estimate": tail estimate of spectral dimension
        - "tail_fraction": tail fraction used
        - "invalid_weight_counts": count of invalid weights discarded
        - "entropy_norm": normalization mode used
        - "legacy_normalization": whether legacy normalization was used
        - "specific_heat_scale": scaling mode for the entropy derivative
    """
    if q <= 0:
        raise ValueError("q must be strictly positive.")
    if steps < 2:
        raise ValueError("steps must be at least 2 to build the entropy profile.")
    if not (0 < tail_fraction <= 1):
        raise ValueError("tail_fraction must be in the interval (0, 1].")

    eigvals = np.asarray(eigenvalues, dtype=typf)
    eps = threshold if threshold is not None else dtype_numerical_precision(typf)
    eigvals = np.where(np.abs(eigvals) > eps, eigvals, typf(0))
    log_N = np.log(typf(num_nodes)) if num_nodes > 1 else typf(1)

    if entropy_norm is None:
        entropy_norm = "standard" if legacy_normalization else "complement"

    tau_grid = np.logspace(t1, t2, steps, dtype=typf)
    renyi_entropy_standard = np.full_like(tau_grid, np.nan, dtype=typf)
    tsallis_entropy_standard = np.full_like(tau_grid, np.nan, dtype=typf)
    log_tau = np.log(tau_grid)

    invalid_counts = 0
    for idx, tau in enumerate(tau_grid):
        if abs(q - 1.0) < 1e-8:
            weights = np.exp(-tau * eigvals)
            Z = np.nansum(weights, dtype=typf)
            if Z == 0:
                continue
            p = weights / Z
            with np.errstate(divide="ignore", invalid="ignore"):
                shannon_raw = -np.nansum(p * np.log(p), dtype=typf)
            shannon_norm = shannon_raw / log_N
            renyi_entropy_standard[idx] = shannon_norm
            tsallis_entropy_standard[idx] = shannon_norm
            continue

        base = 1.0 + (1.0 - q) * tau * eigvals
        if drop_invalid:
            valid = base > 0
            invalid_counts += np.count_nonzero(~valid)
            base = np.where(valid, base, 0.0)
        exponent = typf(1.0 / (q - 1.0))
        with np.errstate(over="ignore", invalid="ignore"):
            eps_i = np.power(base, exponent, dtype=typf)
        Z = np.nansum(eps_i, dtype=typf)
        if Z <= 0:
            continue
        p = eps_i / Z
        mq = np.nansum(np.power(p, q, dtype=typf), dtype=typf)
        if mq <= 0:
            continue
        renyi_raw = (1.0 / (1.0 - q)) * np.log(mq)
        tsallis_raw = (1.0 / (q - 1.0)) * (1.0 - mq)
        renyi_entropy_standard[idx] = renyi_raw / log_N
        tsallis_entropy_standard[idx] = tsallis_raw / log_N

    renyi_entropy = _normalize_entropy_profile(
        renyi_entropy_standard,
        entropy_norm,
    )
    tsallis_entropy = _normalize_entropy_profile(
        tsallis_entropy_standard,
        entropy_norm,
    )

    with np.errstate(invalid="ignore"):
        heat_mode = entropy_norm
        if legacy_normalization and entropy_norm.lower() == "standard":
            heat_mode = "complement"
        heat_source = _normalize_entropy_profile(
            renyi_entropy_standard,
            heat_mode,
        )
        specific_heat = np.gradient(heat_source, log_tau)
        scale = specific_heat_scale.lower()
        if scale == "logn":
            specific_heat = log_N * specific_heat
        elif scale != "none":
            raise ValueError("specific_heat_scale must be 'logN' or 'none'.")

    tail_len = max(1, int(round(tail_fraction * steps)))
    tail_slice = specific_heat[-tail_len:]
    ds_estimate = (
        np.nan if tail_slice.size == 0 else typf(2.0) * np.nanmean(tail_slice, dtype=typf)
    )

    return {
        "q": q,
        "tau": tau_grid,
        "renyi_entropy": renyi_entropy,
        "tsallis_entropy": tsallis_entropy,
        "specific_heat": specific_heat,
        "ds_estimate": ds_estimate,
        "tail_fraction": tail_fraction,
        "invalid_weight_counts": invalid_counts,
        "entropy_norm": entropy_norm,
        "legacy_normalization": legacy_normalization,
        "specific_heat_scale": specific_heat_scale,
    }
