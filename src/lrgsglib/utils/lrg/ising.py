import numpy as np
#
from pathlib import Path
from typing import Sequence, Tuple, List, Optional
#
from numpy.typing import NDArray

try:
    import cupy as cp
    CUPY_AVAILABLE = True
except ImportError:
    CUPY_AVAILABLE = False
#
from ..basic.linalg import compute_recon_ultra, compute_recon
#
def compose_product_ising_states(
    states: Sequence[NDArray[np.int8]]
) -> NDArray[np.int8]:
    """
    Compute the element-wise product of a sequence of Ising spin configurations.

    Each configuration is represented as a NumPy array with entries ±1. All arrays
    must have the same shape.

    Parameters
    ----------
    states : Sequence[NDArray[np.int8]]
        A non-empty sequence of Ising state arrays (dtype `int8`), each entry ±1.

    Returns
    -------
    NDArray[np.int8]
        A new array of the same shape as the inputs, whose entries are the
        product of corresponding entries across all input arrays.

    Raises
    ------
    ValueError
        If `states` is empty or if the input arrays do not all share the same shape.
    """
    if not states:
        raise ValueError("`states` must contain at least one array.")

    first_shape = states[0].shape
    for idx, state in enumerate(states):
        if state.shape != first_shape:
            raise ValueError(
                f"All state arrays must have the same shape: "
                f"array 0 has shape {first_shape}, but array {idx} has shape {state.shape}."
            )

    # Start from a copy of the first state to preserve input data
    result = states[0].copy()
    for state in states[1:]:
        result *= state

    return result

def compose_weighted_ising_state(
    states: Sequence[NDArray[np.int8]],
    weights: Sequence[float]
) -> NDArray[np.int8]:
    """
    Compute a weighted sum of Ising states and project to ±1.

    Parameters
    ----------
    states : Sequence[NDArray[np.int8]]
        Non-empty sequence of arrays with entries ±1, identical shapes.
    weights : Sequence[float]
        Sequence of weights matching the length of `states`.

    Returns
    -------
    NDArray[np.int8]
        Array of same shape as inputs, entries ±1.

    Raises
    ------
    ValueError
        If `states` is empty, lengths of `states` and `weights` differ,
        or input arrays have mismatched shapes.
    """
    if not states:
        raise ValueError("`states` must contain at least one array.")
    if len(states) != len(weights):
        raise ValueError("Length of `weights` must match length of `states`.")

    shape = states[0].shape
    for idx, state in enumerate(states):
        if state.shape != shape:
            raise ValueError(
                f"All state arrays must have the same shape: "
                f"array 0 has shape {shape}, but array {idx} has shape {state.shape}."
            )

    v = sum(w * state for w, state in zip(weights, states))
    sigma = np.sign(v).astype(np.int8)
    sigma[sigma == 0] = 1
    return sigma

def compose_xor_ising_state(
    states: Sequence[NDArray[np.int8]]
) -> NDArray[np.int8]:
    """
    Compute the bitwise XOR of a sequence of Ising spin configurations.

    Each configuration is a NumPy array with entries ±1. We map +1→0, -1→1,
    perform a bitwise XOR across all arrays, then map back 0→+1, 1→−1.

    Parameters
    ----------
    states : Sequence[NDArray[np.int8]]
        Non-empty sequence of Ising state arrays (dtype `int8`), each entry ±1.
        All arrays must share the same shape.

    Returns
    -------
    NDArray[np.int8]
        Array of the same shape as inputs, entries ±1, representing the XOR result.

    Raises
    ------
    ValueError
        If `states` is empty or if the input arrays do not all share the same shape.
    """
    if not states:
        raise ValueError("`states` must contain at least one array.")
    
    base_shape = states[0].shape
    for idx, state in enumerate(states):
        if state.shape != base_shape:
            raise ValueError(
                f"All state arrays must have the same shape: "
                f"array 0 has shape {base_shape}, but array {idx} has shape {state.shape}."
            )

    # Map +1→0, -1→1
    bit_arrays = [(state < 0).astype(np.int8) for state in states]
    # Compute bitwise XOR over the list
    xor_bits = np.bitwise_xor.reduce(bit_arrays)
    # Map back 0→+1, 1→−1
    result = np.where(xor_bits == 0, 1, -1).astype(np.int8)

    return result

def compute_ising_pairwise_energy(
    spins: NDArray[np.int8],
    edges: Sequence[Tuple[int, int, float]],
    use_gpu: bool = False
) -> float:
    """
    Compute the negative weighted sum of pairwise energy for an Ising spin configuration.

    The energy is defined as:
        E = -∑_{(u,v,w) ∈ edges} w * spins[u] * spins[v]

    Parameters
    ----------
    spins : NDArray[np.int8]
        One-dimensional array of spin values, each entry must be +1 or -1.
    edges : Sequence[Tuple[int, int, float]] or Sequence[Tuple[int, int]]
        Sequence of edges, where each edge is either a 2-tuple (node_index_u, node_index_v)
        or a 3-tuple (node_index_u, node_index_v, weight). If 2-tuples are provided,
        or if weight is None in 3-tuples, weights default to 1.0. `node_index_u` and 
        `node_index_v` must be valid indices into `spins`; `weight` is a float or None.
    use_gpu : bool, default False
        If True, use CuPy for GPU acceleration. Requires CuPy to be installed.

    Returns
    -------
    float
        The total energy of the configuration, computed as the negative
        weighted sum of spin–spin products over all edges.

    Raises
    ------
    ValueError
        If `spins` is not one-dimensional, any edge index is out of bounds,
        or CuPy is requested but not available.
    """
    if spins.ndim != 1:
        raise ValueError(f"`spins` must be a 1D array, got shape {spins.shape}")
    
    if use_gpu and not CUPY_AVAILABLE:
        raise ValueError("CuPy is not available. Install CuPy or set use_gpu=False.")

    # Choose the appropriate library
    xp = cp if use_gpu else np

    # Unzip edge components
    if edges:
        # Check if edges have weights (3-tuples) or just node pairs (2-tuples)
        if len(edges[0]) == 3:
            ui, vi, w = zip(*edges)
            # Check if weights are None and replace with 1.0
            w = [1.0 if weight is None else weight for weight in w]
        elif len(edges[0]) == 2:
            ui, vi = zip(*edges)
            w = [1.0] * len(edges)  # Default weights to 1.0
        else:
            raise ValueError("Each edge must be a 2-tuple (u, v) or 3-tuple (u, v, weight)")
    else:
        ui, vi, w = [], [], []

    # Convert to arrays using appropriate library
    ui_arr = xp.array(ui, dtype=int)
    vi_arr = xp.array(vi, dtype=int)
    w_arr = xp.array(w, dtype=float)
    
    # Convert spins to GPU if using CuPy
    spins_xp = xp.asarray(spins) if use_gpu else spins

    n = spins.size
    if ui_arr.size and (ui_arr.min() < 0 or ui_arr.max() >= n):
        raise ValueError("Edge source indices out of bounds")
    if vi_arr.size and (vi_arr.min() < 0 or vi_arr.max() >= n):
        raise ValueError("Edge target indices out of bounds")

    # Compute energy
    energy = -xp.sum(w_arr * spins_xp[ui_arr] * spins_xp[vi_arr])
    
    # Convert back to CPU if using GPU
    if use_gpu:
        energy = float(cp.asnumpy(energy))
    else:
        energy = float(energy)
    
    return energy

def spin_overlap(S1: np.ndarray, S2: np.ndarray, max_overlap: bool = True) -> float:
    """
    Compute the maximum overlap between two spin configurations.

    This function calculates the overlap between two spin matrices by computing
    both the direct overlap and the flipped overlap (accounting for spin inversion).
    The overlap is defined as the average of the product of corresponding spins.

    Parameters
    ----------
    S1 : np.ndarray
        A 2D numpy array representing the first spin configuration.
    S2 : np.ndarray
        A 2D numpy array representing the second spin configuration.

    Returns
    -------
    float
        The maximum overlap value between the two spin configurations.

    Raises
    ------
    AssertionError
        If S1 and S2 do not have the same shape.

    Notes
    -----
    The function assumes that the spin values are numerical (e.g., +1 and -1).
    """
    assert S1.shape == S2.shape, "Matrices must have the same shape!"

    direct_overlap = np.sum(S1 * S2) / S1.size
    if max_overlap:
        flipped_overlap = np.sum(S1 * (-S2)) / S1.size
        return max(direct_overlap, flipped_overlap)
    else:
        return direct_overlap
#
def spin_matching_fraction(S1: np.ndarray, S2: np.ndarray, *, z2: bool = True) -> float:
    """
    Fraction of matching spins between two ±1 configurations.

    Parameters
    ----------
    S1, S2 : np.ndarray
        Arrays with identical shape containing ±1.
    z2 : bool, default True
        If True, maximises over a global spin flip (Z₂ symmetry).

    Returns
    -------
    float
        Matching fraction in [0, 1].
    """
    assert S1.shape == S2.shape, "Arrays must have the same shape."
    same = np.sum(S1 == S2)
    if z2:
        opposite = np.sum(S1 == -S2)
        same = max(same, opposite)
    return same / S1.size
#
def spin_matching_fraction_fromovp(
    S1: np.ndarray, S2: np.ndarray, *, max_overlap: bool = True
) -> float:
    """
    Return the fraction of matching spins between two ±1 configurations,
    using `spin_overlap` for the underlying overlap calculation.

    Parameters
    ----------
    S1, S2 : np.ndarray
        Arrays of identical shape containing ±1 values.
    max_overlap : bool, default True
        Passed to `spin_overlap`; if True, a global Z₂ flip is applied to
        maximise the overlap.

    Returns
    -------
    float
        Matching fraction in [0, 1].  With ``max_overlap=True`` the value
        is ≥ 0.5, since the best orientation is chosen.
    """
    q = spin_overlap(S1, S2, max_overlap=max_overlap)
    return (q + 1.0) / 2.0
#
def compute_spin_overlap_series(spin_vector: np.ndarray, basis: np.ndarray) -> np.ndarray:
    """
    Compute spin overlap between reconstructed and original states.

    Parameters
    -----------
    spin_vector: np.ndarray
        The spin configuration vector.
    basis: np.ndarray
        Array of vectors used for reconstruction.

    Returns
    -----------
    np.ndarray
        Spin overlap series as a function of basis size.
    """
    return np.array([
        spin_overlap(np.sign(compute_recon(spin_vector, basis[:i + 1])), spin_vector)
        for i in range(len(basis)-1)
    ])

def compute_spin_match_series(spin_vector: np.ndarray,
                              basis: np.ndarray,
                              *,
                              z2: bool = True) -> np.ndarray:
    """
    Series of matching fractions between reconstructed and original states.

    Parameters
    ----------
    spin_vector : np.ndarray
        Original spin configuration.
    basis : np.ndarray
        Basis vectors used for reconstruction.
    z2 : bool, default True
        Whether to maximise over a global Z₂ flip.

    Returns
    -------
    np.ndarray
        Matching fraction for each truncation of the basis.
    """
    recon = compute_recon_ultra(spin_vector, basis)
    return np.array([
        spin_matching_fraction(
            np.sign(recon[i]),
            spin_vector,
            z2=z2
        )
        for i in range(len(basis) - 1)
    ])
