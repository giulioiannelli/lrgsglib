/**
 * @file voter_native.cpp
 * @brief Pybind11 wrapper for the voter-model C kernels.
 *
 * In-process voter dynamics on signed CSR graphs -- no file I/O, works with
 * both NX and GT graphs, and is seed-reproducible (unlike the C-subprocess
 * backend, which self-seeds from wall-clock/PID). The actual update reuses the
 * shared voter_model_Nstep kernel, so the Python, C-subprocess and pybind11
 * backends apply identical update logic.
 *
 * Naming convention (Python side):
 *   runlang = "pb_voter"  -> voter_sampling()
 */
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <cstring>
#include <vector>

/* ------------------------------------------------------------------ */
/* Pull in the existing C kernels (reused verbatim).                   */
/* ------------------------------------------------------------------ */
extern "C" {
#include "LRGSG_customs.h"
#include "LRGSG_utils.h"
#include "LRGSG_vm.h"
#include "sfmtrng.h"
}

namespace py = pybind11;

/* ==================================================================
 * Helper: build NodesEdges adjacency from flat numpy CSR arrays.
 * Mirror of GraphCSR in ising_native.cpp.
 *
 *   neigh_indices : int64[total_degree]  (concatenated neighbour lists)
 *   neigh_weights : float64[total_degree](concatenated signed weights)
 *   neigh_ptr     : int64[N+1]           (row pointers into the above)
 * ================================================================== */
struct GraphCSR {
    size_t N;
    std::vector<size_t> nlen;           /* degree of each node */
    std::vector<NodeEdges> node_edges;  /* per-node adjacency  */

    /* Backing storage the NodeEdges point into */
    std::vector<size_t> all_neighbors;
    std::vector<double> all_weights;

    GraphCSR(
        size_t N_,
        const py::array_t<int64_t>& neigh_indices,
        const py::array_t<double>& neigh_weights,
        const py::array_t<int64_t>& neigh_ptr
    ) : N(N_), nlen(N_), node_edges(N_)
    {
        auto ni = neigh_indices.unchecked<1>();
        auto nw = neigh_weights.unchecked<1>();
        auto np_ = neigh_ptr.unchecked<1>();

        size_t total = static_cast<size_t>(ni.size());
        all_neighbors.resize(total);
        all_weights.resize(total);

        for (size_t i = 0; i < total; ++i) {
            all_neighbors[i] = static_cast<size_t>(ni(i));
            all_weights[i] = nw(i);
        }

        for (size_t v = 0; v < N; ++v) {
            size_t start = static_cast<size_t>(np_(v));
            size_t end   = static_cast<size_t>(np_(v + 1));
            nlen[v] = end - start;
            node_edges[v].neighbors = all_neighbors.data() + start;
            node_edges[v].weights   = all_weights.data() + start;
        }
    }
};

/* Seed the global SFMT state directly from seed_val.
 *
 * NOTE: we deliberately do NOT call __set_seed_SFMT() here. That shared helper
 * seeds from a *local* {fixed, fixed, time^pid, ...} array and ignores the
 * global `seed_rand`, so it cannot honour an explicit seed (the C-subprocess
 * backends rely on its time/PID behaviour, which we keep). To get reproducible
 * in-process runs we seed the global `sfmt` (used by RNG_u64) ourselves. */
static void seed_rng(uint64_t seed_val) {
    uint32_t seed_arr[4];
    seed_arr[0] = static_cast<uint32_t>(seed_val & 0xFFFFFFFFu);
    seed_arr[1] = static_cast<uint32_t>((seed_val >> 32) & 0xFFFFFFFFu);
    seed_arr[2] = static_cast<uint32_t>((seed_val * 0x9E3779B97F4A7C15ULL) & 0xFFFFFFFFu);
    seed_arr[3] = static_cast<uint32_t>((seed_val ^ 0xBE11AC1A0ULL) & 0xFFFFFFFFu) | 1u;
    sfmt_init_by_array(&sfmt, seed_arr, 4);
}

/* ==================================================================
 * Voter sampling
 * ==================================================================
 *
 * Runs `n_sweeps` asynchronous voter sweeps (each = N single-node updates
 * via voter_model_Nstep), optionally recording the magnetization before
 * each sweep. Returns (final_spins, magn).
 */
static py::tuple voter_sampling(
    py::array_t<int8_t, py::array::c_style> spins_in,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr,
    size_t n_sweeps,
    uint64_t seed_val,
    bool save_magnetization
) {
    auto s_buf = spins_in.request();
    size_t N = static_cast<size_t>(s_buf.size);

    /* Copy input spins so the caller's array is not mutated */
    std::vector<int8_t> spins(N);
    std::memcpy(spins.data(), s_buf.ptr, N * sizeof(int8_t));

    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    size_t n_rec = save_magnetization ? n_sweeps : 0;
    std::vector<double> magn_out(n_rec);

    {
        py::gil_scoped_release release;

        seed_rng(seed_val);

        spin_tp s = spins.data();
        size_tp nlen_ptr = graph.nlen.data();
        NodesEdges ne = graph.node_edges.data();

        for (size_t step = 0; step < n_sweeps; ++step) {
            if (save_magnetization)
                magn_out[step] = calc_magn(N, s);   /* record before sweep */
            voter_model_Nstep(N, s, nlen_ptr, ne);
        }
    }

    py::array_t<int8_t> out_spins(N);
    std::memcpy(out_spins.mutable_data(), spins.data(), N * sizeof(int8_t));

    py::array_t<double> out_magn(n_rec);
    if (n_rec)
        std::memcpy(out_magn.mutable_data(), magn_out.data(),
                    n_rec * sizeof(double));

    return py::make_tuple(out_spins, out_magn);
}

PYBIND11_MODULE(_voter_native, m) {
    m.doc() = "Native (pybind11) voter-model kernels on signed CSR graphs.";

    m.def("voter_sampling", &voter_sampling,
        R"pbdoc(
Run asynchronous voter dynamics on a signed CSR graph.

Each sweep performs N single-node updates: a node copies sign(w_ij) * s_j of a
uniformly chosen neighbour j (a negative edge weight gives an anti-voter copy).
Isolated (degree-0) nodes keep their state.

Returns
-------
tuple[ndarray, ndarray]
    (final_spins[int8], magnetization_trace[float64]) -- the trace is empty
    when ``save_magnetization`` is False.
        )pbdoc",
        py::arg("spins"), py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("n_sweeps"), py::arg("seed"),
        py::arg("save_magnetization") = true);
}
