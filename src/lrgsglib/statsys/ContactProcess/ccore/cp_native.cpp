/**
 * @file cp_native.cpp
 * @brief Pybind11 wrapper for the SIR contact-process kernel on signed CSR
 *        graphs.
 *
 * In-process SIR contact dynamics -- no file I/O, works with both NX and GT
 * graphs, and is seed-reproducible (unlike the C-subprocess backend, which
 * self-seeds from wall-clock/PID). Reuses the shared cp_infection_rate /
 * cp_recovery_rate helpers from LRGSG_cp.c verbatim, so the rate arithmetic
 * is byte-identical to the ContactProcessSIR executable.
 *
 * The SWEEP SEMANTICS deliberately follow the PYTHON reference loop
 * (ContactProcessBase._sample_py + ContactProcessSIR.ds1step), which is the
 * authoritative path, NOT the C subprocess:
 *
 *   - one sweep visits every node exactly once in a fresh RANDOM PERMUTATION
 *     (np.random.shuffle); the C executable instead draws N nodes uniformly
 *     WITH replacement -- the two differ, and this kernel matches Python;
 *   - the density (mean of the binary state) is recorded BEFORE each sweep
 *     advances the state, so `steps` records cover sweep times 0..steps-1;
 *   - there is no early absorbing-state exit in the Python loop; once the
 *     active count hits zero every rate is 0 (prob = 1 - exp(0) = 0), so the
 *     state is frozen. The kernel skips the dead sweeps as a pure
 *     optimization -- the recorded density trace and final state are
 *     identical to running them.
 *
 * RNG: the global SFMT state, seeded exactly like potts_native.cpp (py<->pb
 * agree only in distribution, not bit-for-bit -- the Python loop draws from
 * numpy).
 *
 * Naming convention (Python side):
 *   runlang = "pb"  -> cp_sir_sampling()
 */
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <vector>

/* ------------------------------------------------------------------ */
/* Pull in the existing C kernels (reused verbatim).                   */
/* ------------------------------------------------------------------ */
extern "C" {
#include "LRGSG_customs.h"
#include "LRGSG_utils.h"
#include "LRGSG_cp.h"
#include "sfmtrng.h"
}

namespace py = pybind11;

/* ==================================================================
 * Helper: build NodeEdges adjacency from flat numpy CSR arrays
 * (neigh_indices/weights concatenated, neigh_ptr row pointers). Mirrors
 * the GraphCSR helper in potts_native.cpp / voter_native.cpp.
 * ================================================================== */
struct GraphCSR {
    size_t N;
    std::vector<size_t> nlen;
    std::vector<NodeEdges> node_edges;
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

/* Seed the global SFMT state directly (reproducible in-process runs).
 * Identical scheme to potts_native.cpp::seed_rng. */
static void seed_rng(uint64_t seed_val) {
    uint32_t seed_arr[4];
    seed_arr[0] = static_cast<uint32_t>(seed_val & 0xFFFFFFFFu);
    seed_arr[1] = static_cast<uint32_t>((seed_val >> 32) & 0xFFFFFFFFu);
    seed_arr[2] = static_cast<uint32_t>((seed_val * 0x9E3779B97F4A7C15ULL) & 0xFFFFFFFFu);
    seed_arr[3] = static_cast<uint32_t>((seed_val ^ 0xBE11AC1A0ULL) & 0xFFFFFFFFu) | 1u;
    sfmt_init_by_array(&sfmt, seed_arr, 4);
}

/* ==================================================================
 * SIR contact-process sampling: `steps` random-permutation sweeps.
 *
 * Per node (binary state, active = 1 / inactive = 0):
 *   active   -> recovers  with prob 1 - exp(-rate),
 *               rate = mu + sum_{w<0, nbr active} (-w)   [cp_recovery_rate]
 *   inactive -> activates with prob 1 - exp(-rate),
 *               rate = sum_{w>0, nbr active} w           [cp_infection_rate]
 *
 * When record_density is set, mean(s) is recorded BEFORE each sweep
 * (matching ContactProcessBase._record's record-then-sweep cadence),
 * giving exactly `steps` records.
 *
 * Returns (final_s[int8], density[float64]); density is empty unless
 * record_density.
 * ================================================================== */
static py::tuple cp_sir_sampling(
    py::array_t<int8_t, py::array::c_style> s_in,
    const py::array_t<int64_t>& nbr_idx,
    const py::array_t<double>&  nbr_w,
    const py::array_t<int64_t>& nbr_ptr,
    size_t N_arg, double mu,
    size_t n_sweeps, uint64_t seed_val,
    bool record_density
) {
    auto s_buf = s_in.request();
    size_t N = static_cast<size_t>(s_buf.size);
    if (N_arg != N)
        throw std::invalid_argument("cp_sir_sampling: N does not match len(s).");

    std::vector<int8_t> state(N);
    std::memcpy(state.data(), s_buf.ptr, N * sizeof(int8_t));

    GraphCSR graph(N, nbr_idx, nbr_w, nbr_ptr);

    std::vector<double> dens_out;
    if (record_density) dens_out.reserve(n_sweeps);

    {
        py::gil_scoped_release release;
        seed_rng(seed_val);

        int8_t* s = state.data();
        size_t sum = 0;
        for (size_t i = 0; i < N; ++i) sum += (s[i] != 0);

        /* Fresh random permutation per sweep (py parity: np.random.shuffle
         * of a persistent node buffer). */
        std::vector<size_t> order(N);
        for (size_t i = 0; i < N; ++i) order[i] = i;

        for (size_t t = 0; t < n_sweeps; ++t) {
            if (record_density)
                dens_out.push_back(static_cast<double>(sum) /
                                   static_cast<double>(N));
            /* Absorbing state (no active site): every rate is 0, the sweep
             * is a guaranteed no-op -- skip it (distribution-identical to
             * the Python loop, which keeps sweeping a frozen state). */
            if (sum == 0) continue;

            for (size_t i = N - 1; i > 0; --i) {
                size_t j = static_cast<size_t>(RNG_u64() % (i + 1));
                size_t tmp = order[i];
                order[i] = order[j];
                order[j] = tmp;
            }

            for (size_t k = 0; k < N; ++k) {
                size_t node = order[k];
                size_t degree = graph.nlen[node];
                NodeEdges edges_node = graph.node_edges[node];
                if (s[node]) {
                    double rate = cp_recovery_rate(node, mu, s, degree,
                                                   edges_node);
                    double prob = 1.0 - exp(-rate);
                    if (prob > 0.0 && RNG_dbl() < prob) {
                        s[node] = 0;
                        --sum;
                    }
                } else {
                    double rate = cp_infection_rate(node, s, degree,
                                                    edges_node);
                    double prob = 1.0 - exp(-rate);
                    if (prob > 0.0 && RNG_dbl() < prob) {
                        s[node] = 1;
                        ++sum;
                    }
                }
            }
        }
    }

    py::array_t<int8_t> out_s(N);
    std::memcpy(out_s.mutable_data(), state.data(), N * sizeof(int8_t));

    py::array_t<double> out_dens(dens_out.size());
    if (!dens_out.empty())
        std::memcpy(out_dens.mutable_data(), dens_out.data(),
                    dens_out.size() * sizeof(double));

    return py::make_tuple(out_s, out_dens);
}

PYBIND11_MODULE(_cp_native, m) {
    m.doc() = "Native (pybind11) SIR contact-process kernel on signed CSR graphs.";

    m.def("cp_sir_sampling", &cp_sir_sampling,
        R"pbdoc(
Run SIR contact-process dynamics on a signed CSR graph.

Parameters
----------
s, nbr_idx, nbr_w, nbr_ptr : ndarray
    Initial int8 binary state (active = 1 / inactive = 0) and the CSR
    adjacency (symmetric, SIGNED weights: positive edges infect, negative
    edges add to the recovery rate).
N, mu : system size and spontaneous recovery rate.
steps, seed : run controls (seed-reproducible).
record_density : record mean(s) once per sweep, BEFORE the sweep.

Returns
-------
tuple[ndarray, ndarray]
    (final_state[int8], density[float64]); the density array has length
    `steps` when record_density, else 0.
        )pbdoc",
        py::arg("s"), py::arg("nbr_idx"), py::arg("nbr_w"),
        py::arg("nbr_ptr"), py::arg("N"), py::arg("mu"),
        py::arg("steps"), py::arg("seed"),
        py::arg("record_density") = false);
}
