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
    "cluster_components",
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

