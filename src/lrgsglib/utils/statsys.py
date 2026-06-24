from collections import Counter

import numpy as np
from numpy.typing import NDArray
#
from .tools import UnionFind
#
__all__ = [
    "boltzmann_factor",
    "find_largest_cluster_circle2D",
    "correlated_binary_sequence_vectorized",
    "edge_sign_arrays",
    "signed_neighbor_arrays",
    "cluster_components",
    "iter_largest_cluster_masks",
    "compute_largest_cluster_masks",
    "cluster_size_distribution",
    "largest_fraction",
]
#
def boltzmann_factor(E: NDArray, T: float, k_B: float = 1.0) -> NDArray:
    """
    Compute the Boltzmann factor for a given array of energies.

    Parameters:
    -----------
    E: NDArray
        Array of energy values.
    T: float
        Temperature of the system (must be positive).
    k_B: float
        Boltzmann constant (default is 1.0).

    Returns:
    -----------
    NDArray
        Array of Boltzmann factors corresponding to the input energies.

    Notes:
    -----------
    The Boltzmann factor is calculated as exp(-E / (k_B * T)).
    """
    from numpy import exp
    if T <= 0:
        raise ValueError("Temperature must be positive.")
    return exp(-E / (k_B * T))

def find_largest_cluster_circle2D(circles, radius):
    """
    Identifies the largest cluster of overlapping circles given their centers and a common radius.

    Parameters:
    -----------
    circles (np.array): A numpy array of tuples where each tuple represents the (x, y) coordinates of a circle's center.
    radius (float): The radius of each circle.

    Returns:
    --------
    list: A list of tuples representing the centers of the circles in the largest cluster.

    Description:
    ------------
    The function utilizes a KDTree for efficient spatial queries to detect overlapping circles.
    It employs a union-find data structure to group and identify clusters of overlapping circles.
    The function returns the largest cluster found.
    """
    from scipy.spatial import KDTree
    from collections import defaultdict
    tree = KDTree(circles)
    uf = UnionFind(len(circles))
    threshold = 2 * radius

    for i in range(len(circles)):
        neighbors = tree.query_ball_point(circles[i], r=threshold)
        for j in neighbors:
            if i != j:
                uf.union(i, j)

    clusters = defaultdict(list)
    for i in range(len(circles)):
        root = uf.find(i)
        clusters[root].append(tuple(circles[i]))

    largest_cluster = max(clusters.values(), key=len)
    return largest_cluster

def correlated_binary_sequence_vectorized(length: int, T: float, J: float = 1.) -> NDArray:
    """
    Generate a random binary sequence with a given flipping probability based on interaction JI and temperature T,
    using a vectorized approach.

    Parameters
    ----------
    length : int
        Length of the sequence.
    J : float
        Interaction strength (coupling constant).
    T : float
        Temperature controlling randomness.

    Returns
    -------
    NDArray
        A binary sequence with flipping probabilities defined by tanh(JI / T).
    """
    from numpy import tanh, random, cumsum, ones
    P_flip = 1-0.5 * (tanh(J / T) + 1)  # Calculate flipping probability
    flips = random.random(size=length) < P_flip
    sequence = ones(length, dtype=int)
    sequence[0] = random.choice([-1, 1])
    flips_cumsum = cumsum(flips)
    sequence = sequence[0] * (-1) ** flips_cumsum
    return sequence

# --------------------------------------------------------------------------- #
# Configuration-domain (cluster) observables                                  #
# --------------------------------------------------------------------------- #
# A *domain* (cluster) is a connected component of the subgraph of **active**
# edges of a spin configuration ``s`` on a (signed) graph. Edge ``(i, j)`` is
# active iff ``b_ij * s_i * s_j > 0``, where the per-edge sign ``b_ij`` selects
# the domain mode:
#   * ``satisfied`` -- ``b_ij = sign(w_ij)``: domains of aligned signed opinion
#     (a negative edge "agrees" when the spins are opposite); reduces to
#     same-spin domains on an unsigned graph.
#   * ``rawspin``   -- ``b_ij = +1``: domains of equal raw spin value.
# These are model-agnostic observables (voter, Ising, Potts, ...): they take the
# ragged per-node neighbour/sign representation (``idx``, ``b``) plus a spin
# vector. The native voter C/pybind kernel mirrors ``cluster_size_distribution``
# exactly and is checked against it for parity (see the voter cluster tests).

def edge_sign_arrays(signs, cluster_mode: str):
    """
    Per-edge ``b_ij`` arrays for the requested ``cluster_mode``.

    Parameters
    ----------
    signs : Sequence[NDArray]
        Ragged per-node array of ``sign(w_ij)`` (the convention returned by
        ``VoterModel._gillespie_neighbors``).
    cluster_mode : str
        ``"satisfied"`` -> edge sign ``sign(w_ij)``; ``"rawspin"`` -> ``+1``
        (edge sign ignored).

    Returns
    -------
    list of NDArray
        Per-node ``int8`` arrays of ``b_ij`` aligned with ``signs``.
    """
    if cluster_mode == "satisfied":
        return [np.asarray(sg, dtype=np.int8) for sg in signs]
    if cluster_mode == "rawspin":
        return [np.ones(len(sg), dtype=np.int8) for sg in signs]
    raise ValueError(f"unknown cluster_mode={cluster_mode!r}")

def signed_neighbor_arrays(sg, cluster_mode: str = "rawspin"):
    """
    Ragged signed-adjacency arrays for cluster analysis, from a signed graph.

    Builds, straight from ``sg``'s own adjacency, the ``(idx, b)`` pair that
    :func:`cluster_components` consumes: ``idx[i]`` holds node ``i``'s neighbour
    indices and ``b[i]`` the matching per-edge signs for ``cluster_mode`` (see
    :func:`edge_sign_arrays`). Engine-agnostic -- needs only ``sg.N`` and
    ``sg.get_neighbors_with_weights`` -- and matches the convention of
    ``VoterModel._gillespie_neighbors`` and the 2D/3D lattice cluster animations.

    Parameters
    ----------
    sg
        A signed graph exposing ``N`` and ``get_neighbors_with_weights(node)``.
    cluster_mode : str
        ``"rawspin"`` or ``"satisfied"`` (see :func:`edge_sign_arrays`).

    Returns
    -------
    (list of NDArray, list of NDArray)
        ``idx`` (``int64`` neighbour indices) and ``b`` (``int8`` edge signs),
        each ragged per node and aligned for :func:`cluster_components`.
    """
    idx, signs = [], []
    for nd in range(sg.N):
        js, sgn = [], []
        for j, w in sg.get_neighbors_with_weights(nd):
            js.append(j)
            sgn.append(-1 if w < 0 else 1)
        idx.append(np.asarray(js, dtype=np.int64))
        signs.append(np.asarray(sgn, dtype=np.int8))
    return idx, edge_sign_arrays(signs, cluster_mode)

def cluster_components(s, idx, b) -> np.ndarray:
    """
    Connected-component label per node via BFS over active edges.

    ``idx[i]`` are node ``i``'s neighbours and ``b[i]`` the matching per-edge
    signs; edge ``(i, idx[i][k])`` is active iff
    ``b[i][k] * s[i] * s[idx[i][k]] > 0``. Components are recomputed from scratch
    in a single ``O(N + E)`` pass -- trivially correct and free of the fragile
    incremental split bookkeeping a maintained partition would need.

    Parameters
    ----------
    s : NDArray
        Spin/opinion vector (``+/-1``).
    idx : Sequence[NDArray]
        Ragged per-node neighbour indices.
    b : Sequence[NDArray]
        Ragged per-node edge signs (see :func:`edge_sign_arrays`).

    Returns
    -------
    NDArray
        ``int64`` array of component labels, one per node.
    """
    n = len(idx)
    label = np.full(n, -1, dtype=np.int64)
    cur = 0
    for start in range(n):
        if label[start] != -1:
            continue
        label[start] = cur
        stack = [start]
        while stack:
            u = stack.pop()
            su = int(s[u])
            ju, bu = idx[u], b[u]
            for k in range(ju.size):
                v = int(ju[k])
                if label[v] == -1 and int(bu[k]) * su * int(s[v]) > 0:
                    label[v] = cur
                    stack.append(v)
        cur += 1
    return label

def flat_signed_edges(idx, b):
    """Flatten ragged ``(idx, b)`` neighbour arrays into deduped *undirected*
    edge arrays ``(I, J, B)`` -- each edge once with ``I < J`` -- for vectorised,
    whole-trajectory clustering (see :func:`iter_largest_cluster_masks`)."""
    sizes = [len(x) for x in idx]
    I = np.repeat(np.arange(len(idx), dtype=np.int64), sizes)
    if I.size:
        J = np.concatenate([np.asarray(x) for x in idx]).astype(np.int64)
        B = np.concatenate([np.asarray(x) for x in b]).astype(np.int8)
    else:
        J = np.empty(0, np.int64)
        B = np.empty(0, np.int8)
    keep = I < J
    return I[keep], J[keep], B[keep]


def iter_largest_cluster_masks(states, idx, b):
    """Yield the boolean mask of the largest active-edge cluster per state.

    Vectorised analogue of ``cluster_components`` + ``argmax``-``bincount``:
    builds the edge list once (:func:`flat_signed_edges`) and labels each frame
    with ``scipy``'s C connected-components (``O(N + E)`` in C, ~15x faster than
    the pure-Python BFS over a 3000-frame trajectory). The active-edge convention
    is identical to :func:`cluster_components` (``b_ij * s_i * s_j > 0``), so the
    *largest-cluster size* matches exactly; the *identity* of the largest cluster
    can differ only under an exact size tie (cosmetic for animation).

    Parameters
    ----------
    states : Iterable[NDArray]
        Trajectory of spin/opinion vectors (``+/-1``).
    idx, b : Sequence[NDArray]
        Ragged neighbour indices and edge signs (see :func:`edge_sign_arrays`).

    Yields
    ------
    NDArray[bool]
        Length-``N`` mask selecting the largest active-edge cluster.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    n = len(idx)
    I, J, B = flat_signed_edges(idx, b)
    ones = np.ones(I.size, dtype=bool)
    for state in states:
        s = np.asarray(state, dtype=np.int8).ravel()
        if I.size:
            act = (B * s[I] * s[J]) > 0
            graph = coo_matrix((ones[act], (I[act], J[act])), shape=(n, n))
            _, label = connected_components(graph, directed=False)
        else:
            label = np.arange(n)
        yield label == np.bincount(label).argmax()


def _largest_cluster_masks_chunk(states, idx, b) -> NDArray:
    """Stack the largest-cluster masks for a block of states (module-level so it
    is picklable for a process pool)."""
    return np.asarray(list(iter_largest_cluster_masks(states, idx, b)), dtype=bool)


def compute_largest_cluster_masks(states, idx, b, *, workers: int | None = None) -> NDArray:
    """Buffer the largest active-edge cluster mask of every frame as one
    ``(n_frames, N)`` boolean array.

    The per-frame connected-components labelling -- not the video encode -- is the
    cost of a largest-cluster movie, so this fans the trajectory out across
    **processes** (true parallelism; the pure-Python ``coo_matrix`` build holds
    the GIL, so threads barely help). Pass the result as ``masks=`` to
    ``lat.animate.largest_cluster`` to skip re-clustering and encode at
    state-movie speed, or reuse it across several renders.

    Parameters
    ----------
    states : Iterable[NDArray]
        Trajectory of spin/opinion vectors (each length ``N`` or ``syshape``).
    idx, b : Sequence[NDArray]
        Ragged neighbour indices and edge signs (see :func:`signed_neighbor_arrays`).
    workers : int, optional
        Process count. ``None`` -> ``min(8, cpu_count)``; ``1`` -> serial (no pool).

    Returns
    -------
    NDArray[bool]
        ``(n_frames, N)`` -- row ``t`` selects the largest cluster of frame ``t``.
    """
    import os

    arr = np.asarray(states)
    arr = arr.reshape(arr.shape[0], -1)
    n = arr.shape[0]
    if n == 0:
        return np.empty((0, len(idx)), dtype=bool)

    w = min(8, os.cpu_count() or 1) if workers is None else max(1, int(workers))
    # Below ~1 chunk's worth, or single worker, the pool overhead isn't worth it.
    if w == 1 or n < 1500:
        return _largest_cluster_masks_chunk(arr, idx, b)

    try:
        from joblib import Parallel, delayed
    except ImportError:
        return _largest_cluster_masks_chunk(arr, idx, b)

    edges = [round(i * n / w) for i in range(w + 1)]
    parts = Parallel(n_jobs=w, backend="loky")(
        delayed(_largest_cluster_masks_chunk)(arr[a:c], idx, b)
        for a, c in zip(edges[:-1], edges[1:]) if c > a
    )
    return np.concatenate(parts, axis=0)


def cluster_size_distribution(s, idx, b) -> Counter:
    """Counter ``{domain_size: number_of_domains}`` for configuration ``s``."""
    label = cluster_components(s, idx, b)
    sizes = np.bincount(label)
    return Counter(int(x) for x in sizes)

def largest_fraction(s_t, idx, b) -> list:
    """
    Largest-domain size / N for each configuration in a trajectory.

    Parameters
    ----------
    s_t : Iterable[NDArray]
        Iterable of spin vectors (e.g. a recorded trajectory ``vm.s_t``).
    idx, b : Sequence[NDArray]
        Neighbour indices and edge signs (see :func:`edge_sign_arrays`).

    Returns
    -------
    list of float
        Giant-domain fraction in ``[0, 1]`` per configuration -- a standard
        coarsening order parameter.
    """
    N = len(idx)
    return [int(np.bincount(cluster_components(np.asarray(s, np.int8), idx, b)).max()) / N
            for s in s_t]

