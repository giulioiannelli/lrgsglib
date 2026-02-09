"""CuPy-based GPU Ising dynamics kernels.

Provides GPU-accelerated Metropolis, Simulated Annealing, and Parallel
Tempering for Ising spin systems on arbitrary graph topologies.

The implementation uses a graph-coloring approach: nodes are partitioned
into independent sets (colors) so that all nodes of the same color can
be updated simultaneously on the GPU without race conditions.
For bipartite graphs (e.g. square lattices), this is the classic
checkerboard decomposition with 2 colors.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

try:
    import cupy as cp
    from cupy import RawKernel

    CUPY_AVAILABLE = True
except ImportError:
    CUPY_AVAILABLE = False
    cp = None  # type: ignore[assignment]

# -----------------------------------------------------------------------
# CUDA kernel: Metropolis update for one color class
# -----------------------------------------------------------------------
_METROPOLIS_KERNEL_SRC = r"""
extern "C" __global__
void metropolis_update(
    signed char* spins,        // [N] spin array (+1/-1)
    const long long* indices,  // CSR neighbor indices
    const double* weights,     // CSR edge weights (signs)
    const long long* ptr,      // CSR row pointers [N+1]
    const int* color_nodes,    // nodes in this color class
    const int n_color,         // size of color class
    const double* randoms,     // [n_color] uniform(0,1)
    const double inv_T         // 1/T (inverse temperature)
) {
    int tid = blockDim.x * blockIdx.x + threadIdx.x;
    if (tid >= n_color) return;

    int node = color_nodes[tid];
    int si = spins[node];

    // Compute local field: sum of J_ij * s_j
    double h_local = 0.0;
    long long start = ptr[node];
    long long end = ptr[node + 1];
    for (long long j = start; j < end; j++) {
        h_local += weights[j] * (double)spins[indices[j]];
    }

    // Energy difference for flipping: dE = 2 * s_i * h_local
    double dE = 2.0 * (double)si * h_local;

    // Metropolis acceptance
    if (dE <= 0.0 || randoms[tid] < exp(-dE * inv_T)) {
        spins[node] = (signed char)(-si);
    }
}
"""

_ENERGY_KERNEL_SRC = r"""
extern "C" __global__
void compute_energy_per_node(
    const signed char* spins,
    const long long* indices,
    const double* weights,
    const long long* ptr,
    double* energy_out,   // [N] per-node energy contribution
    const int N
) {
    int node = blockDim.x * blockIdx.x + threadIdx.x;
    if (node >= N) return;

    double e = 0.0;
    long long start = ptr[node];
    long long end = ptr[node + 1];
    for (long long j = start; j < end; j++) {
        int nb = (int)indices[j];
        if (nb > node) {  // count each edge once
            e -= weights[j] * (double)spins[node] * (double)spins[nb];
        }
    }
    energy_out[node] = e;
}
"""

_MAGNETIZATION_KERNEL_SRC = r"""
extern "C" __global__
void compute_magnetization(
    const signed char* spins,
    double* result,
    const int N
) {
    // Simple reduction - use CuPy's sum for production
    // This kernel is a placeholder; we use cp.sum() instead.
}
"""


def _compile_kernels():
    """Compile CUDA kernels (cached by CuPy)."""
    met_kernel = RawKernel(_METROPOLIS_KERNEL_SRC, "metropolis_update")
    ene_kernel = RawKernel(_ENERGY_KERNEL_SRC, "compute_energy_per_node")
    return met_kernel, ene_kernel


def _greedy_coloring(
    n_nodes: int,
    neigh_indices: np.ndarray,
    neigh_ptr: np.ndarray,
) -> list[np.ndarray]:
    """Greedy graph coloring for parallel Metropolis updates.

    Returns a list of int32 arrays, one per color, containing node indices.
    For bipartite graphs (square lattices) this produces exactly 2 colors.
    """
    colors = np.full(n_nodes, -1, dtype=np.int32)
    n_colors = 0

    for node in range(n_nodes):
        # Find colors used by neighbors
        start, end = int(neigh_ptr[node]), int(neigh_ptr[node + 1])
        nbrs = neigh_indices[start:end]
        used = set()
        for nb in nbrs:
            c = colors[nb]
            if c >= 0:
                used.add(c)
        # Assign smallest available color
        c = 0
        while c in used:
            c += 1
        colors[node] = c
        if c >= n_colors:
            n_colors = c + 1

    # Group nodes by color
    color_classes = []
    for c in range(n_colors):
        nodes = np.where(colors == c)[0].astype(np.int32)
        color_classes.append(nodes)
    return color_classes


def gpu_metropolis_sweep(
    spins_gpu: "cp.ndarray",
    indices_gpu: "cp.ndarray",
    weights_gpu: "cp.ndarray",
    ptr_gpu: "cp.ndarray",
    color_classes_gpu: list["cp.ndarray"],
    inv_T: float,
    met_kernel,
    rng: "cp.random.RandomState",
) -> None:
    """One full Metropolis sweep over all color classes."""
    block = 256
    for cnodes_gpu in color_classes_gpu:
        nc = len(cnodes_gpu)
        if nc == 0:
            continue
        grid = (nc + block - 1) // block
        randoms = rng.rand(nc).astype(cp.float64)
        met_kernel(
            (grid,), (block,),
            (spins_gpu, indices_gpu, weights_gpu, ptr_gpu,
             cnodes_gpu, np.int32(nc), randoms, np.float64(inv_T)),
        )


def gpu_compute_energy(
    spins_gpu: "cp.ndarray",
    indices_gpu: "cp.ndarray",
    weights_gpu: "cp.ndarray",
    ptr_gpu: "cp.ndarray",
    N: int,
    ene_kernel,
) -> float:
    """Compute total Ising energy on GPU."""
    block = 256
    grid = (N + block - 1) // block
    energy_buf = cp.zeros(N, dtype=cp.float64)
    ene_kernel(
        (grid,), (block,),
        (spins_gpu, indices_gpu, weights_gpu, ptr_gpu, energy_buf, np.int32(N)),
    )
    return float(cp.sum(energy_buf))


def gpu_compute_magnetization(spins_gpu: "cp.ndarray", N: int) -> float:
    """Compute magnetization on GPU (raw spin sum, matching Python convention)."""
    return float(cp.sum(spins_gpu).item())


# -----------------------------------------------------------------------
# High-level functions called from IsingDynamics
# -----------------------------------------------------------------------

def cupy_metropolis(
    spins: np.ndarray,
    neigh_indices: np.ndarray,
    neigh_weights: np.ndarray,
    neigh_ptr: np.ndarray,
    T: float,
    n_sweeps: int,
    seed: int,
) -> tuple[np.ndarray, list[float], list[float]]:
    """GPU Metropolis sampling.

    Parameters
    ----------
    spins : ndarray[int8], shape (N,)
    neigh_indices, neigh_weights, neigh_ptr : CSR arrays
    T : float
    n_sweeps : int
    seed : int

    Returns
    -------
    final_spins : ndarray[int8]
    energy_trace : list[float]
    magnetization_trace : list[float]
    """
    if not CUPY_AVAILABLE:
        raise RuntimeError("CuPy not available. Install cupy for GPU backend.")

    N = len(spins)
    inv_T = 1.0 / T if T > 0 else 1e12

    # Graph coloring (CPU, one-time)
    color_classes = _greedy_coloring(N, neigh_indices, neigh_ptr)

    # Compile kernels
    met_kernel, ene_kernel = _compile_kernels()

    # Transfer to GPU
    spins_gpu = cp.asarray(spins, dtype=cp.int8)
    idx_gpu = cp.asarray(neigh_indices, dtype=cp.int64)
    wgt_gpu = cp.asarray(neigh_weights, dtype=cp.float64)
    ptr_gpu = cp.asarray(neigh_ptr, dtype=cp.int64)
    cc_gpu = [cp.asarray(c, dtype=cp.int32) for c in color_classes]

    rng = cp.random.RandomState(seed)

    energy_trace: list[float] = []
    magn_trace: list[float] = []

    for _ in range(n_sweeps):
        gpu_metropolis_sweep(
            spins_gpu, idx_gpu, wgt_gpu, ptr_gpu, cc_gpu,
            inv_T, met_kernel, rng,
        )
        energy_trace.append(
            gpu_compute_energy(spins_gpu, idx_gpu, wgt_gpu, ptr_gpu, N, ene_kernel)
        )
        magn_trace.append(gpu_compute_magnetization(spins_gpu, N))

    final_spins = cp.asnumpy(spins_gpu)
    return final_spins, energy_trace, magn_trace


def cupy_sa(
    spins: np.ndarray,
    neigh_indices: np.ndarray,
    neigh_weights: np.ndarray,
    neigh_ptr: np.ndarray,
    T_init: float,
    T_final: float,
    schedule: str,
    cooling_rate: float,
    steps_per_T: int,
    n_temperatures: int,
    seed: int,
) -> tuple[np.ndarray, list[float], list[float], list[float]]:
    """GPU Simulated Annealing."""
    if not CUPY_AVAILABLE:
        raise RuntimeError("CuPy not available.")

    N = len(spins)

    # Build temperature schedule
    if schedule in ("exponential", "exp"):
        temps = T_init * (cooling_rate ** np.arange(n_temperatures))
    elif schedule == "linear":
        temps = np.linspace(T_init, T_final, n_temperatures)
    else:  # logarithmic
        temps = T_init / (1 + cooling_rate * np.arange(n_temperatures))
    temps = np.clip(temps, T_final, T_init)

    color_classes = _greedy_coloring(N, neigh_indices, neigh_ptr)
    met_kernel, ene_kernel = _compile_kernels()

    spins_gpu = cp.asarray(spins, dtype=cp.int8)
    idx_gpu = cp.asarray(neigh_indices, dtype=cp.int64)
    wgt_gpu = cp.asarray(neigh_weights, dtype=cp.float64)
    ptr_gpu = cp.asarray(neigh_ptr, dtype=cp.int64)
    cc_gpu = [cp.asarray(c, dtype=cp.int32) for c in color_classes]
    rng = cp.random.RandomState(seed)

    energy_trace: list[float] = []
    magn_trace: list[float] = []

    for T in temps:
        inv_T = 1.0 / T if T > 0 else 1e12
        for _ in range(steps_per_T):
            gpu_metropolis_sweep(
                spins_gpu, idx_gpu, wgt_gpu, ptr_gpu, cc_gpu,
                inv_T, met_kernel, rng,
            )
            energy_trace.append(
                gpu_compute_energy(spins_gpu, idx_gpu, wgt_gpu, ptr_gpu, N, ene_kernel)
            )
            magn_trace.append(gpu_compute_magnetization(spins_gpu, N))

    final_spins = cp.asnumpy(spins_gpu)
    return final_spins, energy_trace, magn_trace, temps.tolist()


def cupy_pt(
    spins: np.ndarray,
    neigh_indices: np.ndarray,
    neigh_weights: np.ndarray,
    neigh_ptr: np.ndarray,
    T_ladder: np.ndarray,
    steps_per_exchange: int,
    n_exchanges: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """GPU Parallel Tempering — all replicas updated on GPU simultaneously.

    Returns
    -------
    final_spins : ndarray[int8], shape (N,) — lowest-T replica
    energy_2d : ndarray, shape (n_replicas, n_exchanges)
    magn_2d : ndarray, shape (n_replicas, n_exchanges)
    exchanges : ndarray[bool], shape (n_exchanges, n_replicas-1)
    """
    if not CUPY_AVAILABLE:
        raise RuntimeError("CuPy not available.")

    N = len(spins)
    n_rep = len(T_ladder)

    color_classes = _greedy_coloring(N, neigh_indices, neigh_ptr)
    met_kernel, ene_kernel = _compile_kernels()

    # Transfer graph to GPU (shared across replicas)
    idx_gpu = cp.asarray(neigh_indices, dtype=cp.int64)
    wgt_gpu = cp.asarray(neigh_weights, dtype=cp.float64)
    ptr_gpu = cp.asarray(neigh_ptr, dtype=cp.int64)
    cc_gpu = [cp.asarray(c, dtype=cp.int32) for c in color_classes]

    # Create replicas (all start from same initial state)
    replicas_gpu = [cp.asarray(spins.copy(), dtype=cp.int8) for _ in range(n_rep)]
    rngs = [cp.random.RandomState(seed + i) for i in range(n_rep)]

    energy_2d = np.zeros((n_rep, n_exchanges), dtype=np.float64)
    magn_2d = np.zeros((n_rep, n_exchanges), dtype=np.float64)
    exchanges = np.zeros((n_exchanges, n_rep - 1), dtype=np.bool_)

    cpu_rng = np.random.RandomState(seed + n_rep)

    for ex in range(n_exchanges):
        # Run sweeps for each replica
        for r in range(n_rep):
            inv_T = 1.0 / T_ladder[r] if T_ladder[r] > 0 else 1e12
            for _ in range(steps_per_exchange):
                gpu_metropolis_sweep(
                    replicas_gpu[r], idx_gpu, wgt_gpu, ptr_gpu, cc_gpu,
                    inv_T, met_kernel, rngs[r],
                )

        # Record observables
        for r in range(n_rep):
            energy_2d[r, ex] = gpu_compute_energy(
                replicas_gpu[r], idx_gpu, wgt_gpu, ptr_gpu, N, ene_kernel
            )
            magn_2d[r, ex] = gpu_compute_magnetization(replicas_gpu[r], N)

        # Attempt exchanges (even/odd alternating)
        start = ex % 2
        for i in range(start, n_rep - 1, 2):
            dBeta = (1.0 / T_ladder[i] - 1.0 / T_ladder[i + 1])
            dE = energy_2d[i, ex] - energy_2d[i + 1, ex]
            delta = dBeta * dE
            if delta <= 0 or cpu_rng.random() < np.exp(-delta):
                # Swap replicas (just swap GPU array references)
                replicas_gpu[i], replicas_gpu[i + 1] = (
                    replicas_gpu[i + 1], replicas_gpu[i],
                )
                rngs[i], rngs[i + 1] = rngs[i + 1], rngs[i]
                exchanges[ex, i] = True

    final_spins = cp.asnumpy(replicas_gpu[0])
    return final_spins, energy_2d, magn_2d, exchanges
