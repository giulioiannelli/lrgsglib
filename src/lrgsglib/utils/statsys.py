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
    "cluster_wraps",
    "giant_cluster_spans",
    "consensus_time_stats",
    "survival_curve",
    "interface_density_ensemble",
    "order_parameter_susceptibility",
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


# ---------------------------------------------------------------------------
# Canonical voter / coarsening observables (ensemble-aware)
# ---------------------------------------------------------------------------
def cluster_wraps(mask) -> bool:
    """True iff the ``True`` region of a boolean grid ``mask`` wraps the torus.

    The discrete-percolation *spanning* test: flood-fill the active cells while
    carrying each cell's **unwrapped** integer coordinate; if the fill returns to
    an already-visited cell with a *shifted* unwrapped coordinate along some axis,
    the cluster wraps that axis (it connects to its own periodic image). This is
    the standard wrapping criterion -- strictly stronger than a naive
    "touches both opposite faces", which can false-positive.

    ``mask`` carries one axis per lattice dimension (e.g. a per-node giant-cluster
    mask reshaped to ``sg.syshape``); neighbours are the von Neumann ``+/-1`` per
    axis, so this is exact for square / cubic lattices. N-D generalisation of the
    2D test used in the ``ipy/labs/VM`` giant-cluster lab.

    Parameters
    ----------
    mask : NDArray[bool]
        Boolean array, one axis per lattice dimension.

    Returns
    -------
    bool
        ``True`` iff the masked region wraps along at least one axis.
    """
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        return False
    shape = mask.shape
    ndim = mask.ndim
    visited = np.zeros(shape, dtype=bool)
    # unwrapped[d] holds the un-modded coordinate along axis d for each cell.
    unwrapped = np.zeros((ndim,) + shape, dtype=np.int64)
    wrap = [False] * ndim
    offsets = []
    for d in range(ndim):
        for step in (1, -1):
            off = [0] * ndim
            off[d] = step
            offsets.append(tuple(off))
    for seed in np.argwhere(mask):
        seed = tuple(int(x) for x in seed)
        if visited[seed]:
            continue
        visited[seed] = True
        for d in range(ndim):
            unwrapped[(d,) + seed] = seed[d]
        stack = [seed]
        while stack:
            cur = stack.pop()
            cur_uw = [int(unwrapped[(d,) + cur]) for d in range(ndim)]
            for off in offsets:
                nb = tuple((cur[d] + off[d]) % shape[d] for d in range(ndim))
                if not mask[nb]:
                    continue
                nb_uw = [cur_uw[d] + off[d] for d in range(ndim)]
                if not visited[nb]:
                    visited[nb] = True
                    for d in range(ndim):
                        unwrapped[(d,) + nb] = nb_uw[d]
                    stack.append(nb)
                else:
                    for d in range(ndim):
                        if nb_uw[d] != int(unwrapped[(d,) + nb]):
                            wrap[d] = True
        if all(wrap):
            break
    return any(wrap)


def giant_cluster_spans(s, idx, b, shape) -> bool:
    """True iff the largest active-edge cluster of ``s`` spans (wraps) the torus.

    Combines the engine-agnostic clustering (:func:`cluster_components`, exact for
    any geometry) with the grid wrapping test (:func:`cluster_wraps`): the giant
    component is identified from the real adjacency ``(idx, b)``, then its boolean
    mask is reshaped to ``shape`` and tested for torus wrapping. The reshape
    assumes row-major grid-ordered nodes and a square / cubic lattice (the
    convention of ``LatticeND`` / ``sg.syshape``).

    Parameters
    ----------
    s : NDArray
        Spin / opinion vector (``+/-1``).
    idx, b : Sequence[NDArray]
        Ragged neighbour indices and edge signs (see :func:`signed_neighbor_arrays`).
    shape : tuple of int
        Lattice grid shape (e.g. ``sg.syshape``); ``prod(shape) == len(idx)``.

    Returns
    -------
    bool
        ``True`` iff the giant cluster wraps the torus along some axis.
    """
    label = cluster_components(np.asarray(s, np.int8), idx, b)
    if label.size == 0:
        return False
    giant = int(np.bincount(label).argmax())
    mask = (label == giant).reshape(shape)
    return cluster_wraps(mask)


def consensus_time_stats(absorbed_times, n_steps) -> dict:
    """Exit / consensus-time summary over an ensemble of runs.

    ``absorbed_times`` is one entry per run: the sweep index at which the run
    reached an absorbing (consensus / frozen) state, or ``None`` if it never did
    (right-censored at ``n_steps``). Statistics (mean / std / median) are taken
    over the runs that *did* absorb; ``p_consensus`` is the absorbed fraction and
    ``censored`` counts the rest.

    Parameters
    ----------
    absorbed_times : Sequence[int | None]
        Per-run absorption sweep (``None`` = never absorbed).
    n_steps : int
        Run length (the censoring horizon).

    Returns
    -------
    dict
        ``{n, p_consensus, censored, mean, std, median, n_steps}``; the time
        statistics are ``nan`` when no run absorbed.
    """
    times = list(absorbed_times)
    n = len(times)
    absorbed = np.asarray([int(t) for t in times if t is not None], dtype=float)
    n_abs = absorbed.size
    nan = float("nan")
    return {
        "n": n,
        "p_consensus": (n_abs / n) if n else nan,
        "censored": n - n_abs,
        "mean": float(absorbed.mean()) if n_abs else nan,
        "std": float(absorbed.std(ddof=0)) if n_abs else nan,
        "median": float(np.median(absorbed)) if n_abs else nan,
        "n_steps": int(n_steps),
    }


def survival_curve(absorbed_times, n_steps) -> NDArray:
    """Survival function ``S(t)`` = fraction of runs not yet absorbed by sweep ``t``.

    A run that absorbs at sweep ``a`` is alive for ``t < a`` and dead afterwards;
    a run that never absorbs (``None``) is alive at every ``t``. The result is
    monotone non-increasing with ``S(0) = 1`` (no run absorbs before the first
    sweep).

    Parameters
    ----------
    absorbed_times : Sequence[int | None]
        Per-run absorption sweep (``None`` = never absorbed).
    n_steps : int
        Length of the returned curve.

    Returns
    -------
    NDArray
        ``float`` array of length ``n_steps`` with ``S(t) in [0, 1]``.
    """
    n_steps = int(n_steps)
    times = list(absorbed_times)
    n = len(times)
    if n == 0:
        return np.ones(n_steps, dtype=float)
    t = np.arange(n_steps)
    alive = np.zeros(n_steps, dtype=float)
    for a in times:
        if a is None:
            alive += 1.0
        else:
            alive += (t < int(a)).astype(float)
    return alive / n


def interface_density_ensemble(rho_series) -> tuple:
    """Ensemble mean and std of the interface-density coarsening curve ``rho(t)``.

    Runs that froze early give shorter ``rho(t)`` series; each is **forward-filled**
    at its last value (an absorbed run keeps its terminal ``rho``, typically 0 at
    consensus) up to the longest run, then the ensemble mean / std are taken per
    sweep.

    Parameters
    ----------
    rho_series : Sequence[Sequence[float]]
        One interface-density series per run (e.g. each
        ``VoterModel.interface_density_series()``).

    Returns
    -------
    (NDArray, NDArray)
        ``(mean_rho_t, std_rho_t)``, both length = the longest run (empty if no
        non-empty series were given).
    """
    series = [np.asarray(r, dtype=float) for r in rho_series if len(r)]
    if not series:
        empty = np.empty(0, dtype=float)
        return empty, empty
    max_len = max(s.size for s in series)
    stacked = np.empty((len(series), max_len), dtype=float)
    for i, s in enumerate(series):
        if s.size < max_len:
            s = np.concatenate([s, np.full(max_len - s.size, s[-1])])
        stacked[i] = s
    return stacked.mean(axis=0), stacked.std(axis=0, ddof=0)


def order_parameter_susceptibility(op_samples, n_sites=None) -> tuple:
    """Ensemble mean and susceptibility of an order parameter across disorder.

    ``chi = n_sites * Var(op)`` (population variance ``<op^2> - <op>^2``), the
    standard finite-system fluctuation susceptibility. Pass ``n_sites`` (the
    system size ``N``) for the finite-size-scaling normalisation; with
    ``n_sites=None`` the bare variance is returned (``chi = Var(op)``).

    Parameters
    ----------
    op_samples : Sequence[float]
        Order-parameter value per disorder realisation.
    n_sites : int, optional
        System size ``N`` for the ``N * Var`` normalisation; default ``None``.

    Returns
    -------
    (float, float)
        ``(mean, chi)``; both ``nan`` for an empty sample.
    """
    arr = np.asarray(list(op_samples), dtype=float)
    if arr.size == 0:
        return float("nan"), float("nan")
    var = float(arr.var(ddof=0))
    chi = var * (1.0 if n_sites is None else float(n_sites))
    return float(arr.mean()), chi

