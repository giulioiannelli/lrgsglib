"""
Quantum Propagator Extensions for Signed Laplacian RG

This module implements quantum propagators e^(-itL) for signed graphs,
extending the classical propagators e^(-τL) with quantum interference effects.

Mathematical Background
-----------------------
Classical propagator: U_C(τ) = exp(-τL) - diffusive, real symmetric
Quantum propagator: U_Q(t) = exp(-itL) - coherent, complex unitary

Key differences:
- Classical: eigenvalues decay as exp(-τλ) → spreads diffusively
- Quantum: eigenvalues oscillate as exp(-itλ) → coherent interference

Quantum density matrix evolution:
    ρ(t) = U(t) ρ₀ U†(t)

Quantum probability distribution (localized initial state at node i):
    P_j(t) = |⟨j|U(t)|i⟩|² = |∑ₙ exp(-itλₙ) vₙ(j) vₙ(i)|²

Interference effects arise from cross-terms between eigenmodes:
    P(j,t) = ∑ₙ |vₙ(j)|²|vₙ(i)|²  +  2∑ₙ<ₘ cos((λₙ-λₘ)t) vₙ(j)vₙ(i)vₘ(j)vₘ(i)
             ^^^^^^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
             Classical (incoherent)     Quantum interference (coherent)
"""

import numpy as np
from typing import Tuple, Optional, Callable, Union, Dict
from numpy.typing import NDArray
from scipy.linalg import expm
from scipy.sparse import csr_matrix, csr_array
from scipy.sparse.linalg import expm as sparse_expm

from ..basic import dtype_numerical_precision

__all__ = [
    "compute_quantum_propagator_spectral",
    "compute_quantum_propagator_matrix",
    "quantum_density_matrix_evolution",
    "quantum_probability_distribution",
    "compute_quantum_coherence",
    "von_neumann_entropy",
    "quantum_classical_divergence",
    "compute_interference_visibility",
    "compute_quantum_observables_from_eigenvalues",
]


def compute_quantum_propagator_spectral(
    eigenvalues: NDArray,
    eigenvectors: NDArray,
    t: float,
    full_matrix: bool = False
) -> Union[Tuple[NDArray, NDArray], NDArray]:
    """
    Compute quantum propagator U(t) = exp(-i*t*L) using spectral decomposition.

    Parameters
    ----------
    eigenvalues : NDArray, shape (N,)
        Eigenvalues of signed Laplacian (sorted ascending).
    eigenvectors : NDArray, shape (N, N)
        Eigenvectors of signed Laplacian (column-major: each column is an
        eigenvector).
    t : float
        Quantum time parameter (real-valued).
    full_matrix : bool, optional
        If True, return full U(t) matrix. If False, return diagonal phases
        and eigenvectors for efficient computation (default: False).

    Returns
    -------
    If full_matrix=False (default):
        phases : NDArray, shape (N,), dtype=complex128
            Diagonal elements: exp(-i*t*λ)
        eigenvectors : NDArray, shape (N, N)
            Eigenvectors (for external reconstruction)
    If full_matrix=True:
        U_t : NDArray, shape (N, N), dtype=complex128
            Quantum propagator matrix U(t) = V @ diag(exp(-i*t*λ)) @ V^†

    Notes
    -----
    The quantum propagator is unitary: U(t)^† U(t) = I

    Mathematical formula::

        U(t) = ∑ₙ exp(-i*t*λₙ) \\|vₙ⟩⟨vₙ\\|
             = V · diag(exp(-i*t*eigenvalues)) · V^†

    Examples
    --------

    >>> import numpy as np  # doctest: +SKIP
    >>> from lrgsglib import Lattice2D  # doctest: +SKIP
    >>> l = Lattice2D(side=10, pflip=0.1)  # doctest: +SKIP
    >>> l.flip_random_fract_edges()  # doctest: +SKIP
    >>> l.compute_k_eigvV()  # doctest: +SKIP
    >>> U_t = compute_quantum_propagator_spectral(  # doctest: +SKIP
    ...     l.eigv, l.eigV, t=1.0, full_matrix=True
    ... )  # doctest: +SKIP
    >>> identity = U_t @ U_t.conj().T  # doctest: +SKIP
    >>> print(np.allclose(identity, np.eye(l.N)))  # doctest: +SKIP
    True
    """
    N = len(eigenvalues)

    # Compute quantum phases: exp(-i*t*λ)
    phases = np.exp(-1j * t * eigenvalues)

    if not full_matrix:
        return phases, eigenvectors

    # Reconstruct full propagator: U(t) = V @ diag(phases) @ V^†
    # Efficient computation: U = V @ (phases[:, None] * V.T)
    U_t = eigenvectors @ (phases[:, None] * eigenvectors.T)

    return U_t


def compute_quantum_propagator_matrix(
    L: Union[NDArray, csr_matrix],
    t: float,
) -> NDArray:
    """
    Compute quantum propagator U(t) = exp(-i*t*L) using direct matrix
    exponential.

    Parameters
    ----------
    L : NDArray or csr_matrix
        Signed Laplacian matrix (dense or sparse).
    t : float
        Quantum time parameter.

    Returns
    -------
    U_t : NDArray, dtype=complex128
        Quantum propagator matrix.

    Notes
    -----
    This function uses scipy.linalg.expm or scipy.sparse.linalg.expm.
    For repeated evaluations, use compute_quantum_propagator_spectral instead
    (precompute eigendecomposition once, then vary t efficiently).

    Examples
    --------

    >>> from lrgsglib import Lattice2D  # doctest: +SKIP
    >>> l = Lattice2D(side=10)  # doctest: +SKIP
    >>> L = l.signed_laplacian_matrix(l.G)  # doctest: +SKIP
    >>> U_t = compute_quantum_propagator_matrix(L.toarray(), t=1.0)  # doctest: +SKIP
    >>> print(np.allclose(U_t @ U_t.conj().T, np.eye(l.N)))  # doctest: +SKIP
    True
    """
    if isinstance(L, (csr_matrix, csr_array)):
        # Sparse computation
        U_t = sparse_expm(-1j * t * L)
        return U_t.toarray()
    else:
        # Dense computation
        U_t = expm(-1j * t * L)
        return U_t


def quantum_density_matrix_evolution(
    eigenvalues: NDArray,
    eigenvectors: NDArray,
    t: float,
    rho0: Optional[NDArray] = None,
    init_type: str = "uniform",
    init_node: Optional[int] = None,
) -> NDArray:
    """
    Evolve quantum density matrix: ρ(t) = U(t) ρ₀ U†(t).

    Parameters
    ----------
    eigenvalues : NDArray
        Laplacian eigenvalues.
    eigenvectors : NDArray
        Laplacian eigenvectors (column-major).
    t : float
        Quantum time.
    rho0 : NDArray, optional
        Initial density matrix. If None, generated from init_type.
    init_type : str, optional
        Type of initial state: "uniform", "localized", "thermal"
        (default: "uniform").
    init_node : int, optional
        Node index for localized initial state (required if
        init_type="localized").

    Returns
    -------
    rho_t : NDArray, shape (N, N), dtype=complex128
        Evolved density matrix at time t.

    Notes
    -----
    The density matrix is Hermitian: ρ = ρ†
    Trace is preserved: Tr(ρ(t)) = 1

    Examples
    --------

    >>> from lrgsglib import Lattice2D  # doctest: +SKIP
    >>> l = Lattice2D(side=10)  # doctest: +SKIP
    >>> l.compute_k_eigvV()  # doctest: +SKIP
    >>> rho_t = quantum_density_matrix_evolution(  # doctest: +SKIP
    ...     l.eigv, l.eigV, t=1.0,
    ...     init_type="localized", init_node=50
    ... )  # doctest: +SKIP
    >>> print(np.isclose(np.trace(rho_t), 1.0))  # doctest: +SKIP
    True
    >>> print(np.allclose(rho_t, rho_t.conj().T))  # doctest: +SKIP
    True
    """
    N = len(eigenvalues)

    # Generate initial density matrix if not provided
    if rho0 is None:
        if init_type == "uniform":
            rho0 = np.eye(N, dtype=complex) / N
        elif init_type == "localized":
            if init_node is None:
                raise ValueError(
                    "init_node required for localized initial state"
                )
            rho0 = np.zeros((N, N), dtype=complex)
            rho0[init_node, init_node] = 1.0
        elif init_type == "thermal":
            # Use classical thermal state as initial condition
            tau_init = 1.0 / np.max(eigenvalues) if np.max(eigenvalues) > 0 else 1.0
            thermal_probs = np.exp(-tau_init * eigenvalues)
            thermal_probs /= np.sum(thermal_probs)
            rho0 = eigenvectors @ np.diag(thermal_probs) @ eigenvectors.T
        else:
            raise ValueError(f"Unknown init_type: {init_type}")

    # Compute propagator phases
    phases, V = compute_quantum_propagator_spectral(
        eigenvalues, eigenvectors, t, full_matrix=False
    )

    # Evolve density matrix: ρ(t) = U(t) ρ₀ U†(t)
    # Efficient computation using eigendecomposition
    U_t = V @ (phases[:, None] * V.T)
    rho_t = U_t @ rho0 @ U_t.conj().T

    return rho_t


def quantum_probability_distribution(
    eigenvalues: NDArray,
    eigenvectors: NDArray,
    t: float,
    init_node: int,
) -> NDArray:
    """
    Compute quantum walk probability distribution P_j(t) = |⟨j|U(t)|i⟩|².

    Parameters
    ----------
    eigenvalues : NDArray
        Laplacian eigenvalues.
    eigenvectors : NDArray
        Laplacian eigenvectors.
    t : float
        Quantum time.
    init_node : int
        Initial node position.

    Returns
    -------
    probabilities : NDArray, shape (N,)
        Quantum probability at each node: P_j(t) = |amplitude_j(t)|²

    Notes
    -----
    The quantum amplitude at node j starting from node i:
        ψ_j(t) = ⟨j|U(t)|i⟩ = ∑ₙ exp(-i*t*λₙ) vₙ(j) vₙ(i)

    Probability: P_j(t) = |ψ_j(t)|²

    Probability is conserved: ∑_j P_j(t) = 1

    Examples
    --------
    >>> from lrgsglib import Lattice2D
    >>> l = Lattice2D(side=10)
    >>> l.compute_k_eigvV()
    >>> P_t = quantum_probability_distribution(
    ...     l.eigv, l.eigV, t=1.0, init_node=50
    ... )
    >>> print(np.isclose(np.sum(P_t), 1.0))
    True
    """
    N = len(eigenvalues)

    # Compute quantum phases
    phases = np.exp(-1j * t * eigenvalues)

    # Compute quantum amplitude: ψ_j(t) = ∑ₙ exp(-i*t*λₙ) vₙ(j) vₙ(i)
    # eigenvectors[:, n] is the n-th eigenvector
    amplitudes = np.zeros(N, dtype=complex)
    for n in range(N):
        amplitudes += phases[n] * eigenvectors[:, n] * eigenvectors[init_node, n]

    # Probability is modulus squared
    probabilities = np.abs(amplitudes) ** 2

    return probabilities


def compute_quantum_coherence(rho_t: NDArray) -> float:
    """
    Compute quantum coherence as sum of off-diagonal magnitudes.

    Parameters
    ----------
    rho_t : NDArray
        Quantum density matrix (complex Hermitian).

    Returns
    -------
    coherence : float
        Total quantum coherence: C(ρ) = ∑_{i≠j} abs(ρ_ij)

    Notes
    -----
    Coherence measures the "quantumness" of the state:
    - C = 0: Classical (diagonal density matrix)
    - C > 0: Quantum (off-diagonal coherences present)

    For a pure state |ψ⟩ = ∑_i c_i |i⟩:
    - Maximally coherent: C ≈ N-1 (equal superposition)
    - Minimally coherent: C = 0 (localized state after long time)

    Examples
    --------
    >>> import numpy as np
    >>> # Classical state (diagonal)
    >>> rho_classical = np.diag([0.5, 0.5])
    >>> print(compute_quantum_coherence(rho_classical))
    0.0
    >>> # Quantum superposition
    >>> rho_superposition = np.array([[0.5, 0.5], [0.5, 0.5]])
    >>> print(compute_quantum_coherence(rho_superposition))
    1.0
    """
    N = rho_t.shape[0]
    coherence = 0.0
    for i in range(N):
        for j in range(N):
            if i != j:
                coherence += np.abs(rho_t[i, j])
    return coherence


def von_neumann_entropy(
    rho_t: NDArray,
    numerical_threshold: Optional[float] = None
) -> float:
    """
    Compute von Neumann entropy: S_vN(ρ) = -Tr(ρ log ρ).

    Parameters
    ----------
    rho_t : NDArray
        Quantum density matrix.
    numerical_threshold : float, optional
        Eigenvalues below this threshold treated as zero. If None, uses
        machine precision for complex128 (default: None).

    Returns
    -------
    entropy : float
        von Neumann entropy in nats (use log base 2 for bits).

    Notes
    -----
    For pure states (\\|ψ⟩⟨ψ\\|): S_vN = 0
    For maximally mixed state (I/N): S_vN = log(N)
    General bounds: 0 ≤ S_vN ≤ log(N)

    The von Neumann entropy quantifies the mixedness (classical uncertainty)
    of a quantum state. For bipartite systems, it also measures entanglement.

    Examples
    --------
    >>> import numpy as np
    >>> # Pure state
    >>> rho_pure = np.array([[1.0, 0.0], [0.0, 0.0]], dtype=complex)
    >>> print(von_neumann_entropy(rho_pure))
    0.0
    >>> # Maximally mixed state
    >>> rho_mixed = np.eye(2, dtype=complex) / 2
    >>> print(np.isclose(von_neumann_entropy(rho_mixed), np.log(2)))
    True
    """
    if numerical_threshold is None:
        numerical_threshold = dtype_numerical_precision(np.complex128)

    # Diagonalize density matrix
    eigenvals = np.linalg.eigvalsh(rho_t)

    # Filter small eigenvalues
    eigenvals = eigenvals[eigenvals > numerical_threshold]

    # Compute entropy: -∑ λ log(λ)
    with np.errstate(divide='ignore', invalid='ignore'):
        entropy = -np.sum(eigenvals * np.log(eigenvals))

    return entropy


def quantum_classical_divergence(
    eigenvalues: NDArray,
    eigenvectors: NDArray,
    t: float,
    tau: Optional[float] = None,
    init_node: int = 0,
    metric: str = "frobenius",
) -> float:
    """
    Measure divergence between quantum and classical evolution.

    Parameters
    ----------
    eigenvalues : NDArray
        Laplacian eigenvalues.
    eigenvectors : NDArray
        Laplacian eigenvectors.
    t : float
        Quantum time.
    tau : float, optional
        Classical diffusion time (defaults to t for comparison).
    init_node : int, optional
        Initial node (default: 0).
    metric : str, optional
        Distance metric: "frobenius", "trace", "max" (default: "frobenius").

    Returns
    -------
    divergence : float
        Distance between quantum and classical density matrices.

    Notes
    -----
    Computes D = ||ρ_Q(t) - ρ_C(τ)||

    When divergence is large, quantum interference effects dominate.
    At long times, quantum and classical may converge (depends on graph).

    The classical density matrix here uses a thermal state:
        ρ_C = exp(-τL) / Tr(exp(-τL))

    Examples
    --------
    >>> from lrgsglib import Lattice2D
    >>> l = Lattice2D(side=10)
    >>> l.compute_k_eigvV()
    >>> div = quantum_classical_divergence(
    ...     l.eigv, l.eigV, t=1.0, init_node=50
    ... )
    >>> print(div > 0)
    True
    """
    if tau is None:
        tau = t

    # Compute quantum density matrix
    rho_quantum = quantum_density_matrix_evolution(
        eigenvalues, eigenvectors, t, init_type="localized",
        init_node=init_node
    )

    # Compute classical density matrix (thermal)
    # Avoid overflow/underflow
    max_eigenval = np.max(eigenvalues)
    if max_eigenval > 0:
        shifted_eigenvals = eigenvalues - np.min(eigenvalues)
        classical_probs = np.exp(-tau * shifted_eigenvals)
    else:
        classical_probs = np.exp(-tau * eigenvalues)

    classical_probs /= np.sum(classical_probs)
    rho_classical = eigenvectors @ np.diag(classical_probs) @ eigenvectors.T

    # Compute difference
    diff = rho_quantum - rho_classical

    if metric == "frobenius":
        divergence = np.linalg.norm(diff, 'fro')
    elif metric == "trace":
        divergence = np.trace(np.abs(diff)).real
    elif metric == "max":
        divergence = np.max(np.abs(diff))
    else:
        raise ValueError(f"Unknown metric: {metric}")

    return divergence


def compute_interference_visibility(
    eigenvalues: NDArray,
    eigenvectors: NDArray,
    t_array: NDArray,
    node_i: int,
    node_j: int,
) -> Tuple[NDArray, float]:
    """
    Compute visibility of quantum interference fringes between two nodes.

    Parameters
    ----------
    eigenvalues : NDArray
        Laplacian eigenvalues.
    eigenvectors : NDArray
        Laplacian eigenvectors.
    t_array : NDArray
        Array of time values to scan.
    node_i : int
        First node (source).
    node_j : int
        Second node (target).

    Returns
    -------
    transition_prob : NDArray
        |⟨j|U(t)|i⟩|² for each t in t_array.
    visibility : float
        Fringe visibility: V = (P_max - P_min) / (P_max + P_min)

    Notes
    -----
    High visibility indicates strong quantum interference.
    Visibility ranges from 0 (no interference) to 1 (perfect interference).

    The transition probability oscillates due to beating between eigenfrequencies:
        ω_nm = λ_n - λ_m

    Examples
    --------
    >>> from lrgsglib import Lattice2D
    >>> import numpy as np
    >>> l = Lattice2D(side=10)
    >>> l.compute_k_eigvV()
    >>> t_vals = np.linspace(0, 10, 100)
    >>> prob_t, vis = compute_interference_visibility(
    ...     l.eigv, l.eigV, t_vals, node_i=50, node_j=60
    ... )
    >>> print(0 <= vis <= 1)
    True
    """
    transition_prob = np.zeros(len(t_array))

    for idx, t in enumerate(t_array):
        # Compute quantum amplitude: ⟨j|U(t)|i⟩
        phases = np.exp(-1j * t * eigenvalues)
        amplitude = np.sum(
            phases * eigenvectors[node_j, :] * eigenvectors[node_i, :]
        )
        transition_prob[idx] = np.abs(amplitude) ** 2

    # Compute visibility
    P_max = np.max(transition_prob)
    P_min = np.min(transition_prob)

    if P_max + P_min == 0:
        visibility = 0.0
    else:
        visibility = (P_max - P_min) / (P_max + P_min)

    return transition_prob, visibility


def compute_quantum_observables_from_eigenvalues(
    eigenvalues: NDArray,
    eigenvectors: NDArray,
    num_nodes: int,
    t_array: NDArray,
    init_type: str = "localized",
    init_node: int = 0,
) -> Dict[str, NDArray]:
    """
    Compute quantum observables over time array (batch computation).

    Parameters
    ----------
    eigenvalues : NDArray
        Laplacian eigenvalues.
    eigenvectors : NDArray
        Laplacian eigenvectors.
    num_nodes : int
        Number of nodes in graph.
    t_array : NDArray
        Time values to compute observables.
    init_type : str, optional
        Initial state type (default: "localized").
    init_node : int, optional
        Initial node for localized state (default: 0).

    Returns
    -------
    observables : dict
        Dictionary containing time-series of:
        - "time": Time array (copy of t_array)
        - "von_neumann_entropy": von Neumann entropy S_vN(t)
        - "coherence": Total quantum coherence C(t)
        - "purity": Purity Tr(ρ²) (1 for pure, 1/N for maximally mixed)
        - "participation_ratio": Inverse participation ratio
        - "probability_distribution": P_j(t) for all nodes j (shape: n_times × N)

    Notes
    -----
    This function is optimized for batch computation over many time points.
    All observables are computed in a single pass through the time array.

    Examples
    --------
    >>> from lrgsglib import Lattice2D
    >>> import numpy as np
    >>> l = Lattice2D(side=10)
    >>> l.compute_k_eigvV()
    >>> t_vals = np.linspace(0, 10, 100)
    >>> obs = compute_quantum_observables_from_eigenvalues(
    ...     l.eigv, l.eigV, l.N, t_vals, init_node=50
    ... )
    >>> print(obs.keys())
    dict_keys(['time', 'von_neumann_entropy', 'coherence', 'purity',
               'participation_ratio', 'probability_distribution'])
    """
    N = num_nodes
    n_times = len(t_array)

    # Initialize storage
    entropy_profile = np.zeros(n_times)
    coherence_profile = np.zeros(n_times)
    purity_profile = np.zeros(n_times)
    participation_ratio = np.zeros(n_times)
    prob_distributions = np.zeros((n_times, N))

    for idx, t in enumerate(t_array):
        # Evolve density matrix
        rho_t = quantum_density_matrix_evolution(
            eigenvalues, eigenvectors, t,
            init_type=init_type, init_node=init_node
        )

        # Compute observables
        entropy_profile[idx] = von_neumann_entropy(rho_t)
        coherence_profile[idx] = compute_quantum_coherence(rho_t)

        # Purity: Tr(ρ²)
        purity_profile[idx] = np.trace(rho_t @ rho_t).real

        # Probability distribution (diagonal elements)
        prob_distributions[idx, :] = np.abs(np.diag(rho_t))

        # Participation ratio: 1 / ∑_i P_i²
        probs = prob_distributions[idx, :]
        participation_ratio[idx] = 1.0 / np.sum(probs ** 2)

    observables = {
        "time": t_array.copy(),
        "von_neumann_entropy": entropy_profile,
        "coherence": coherence_profile,
        "purity": purity_profile,
        "participation_ratio": participation_ratio,
        "probability_distribution": prob_distributions,
    }

    return observables
