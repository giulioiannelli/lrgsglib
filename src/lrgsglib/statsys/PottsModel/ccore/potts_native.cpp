/**
 * @file potts_native.cpp
 * @brief Pybind11 wrapper for the q-state Potts Metropolis C kernel.
 *
 * In-process Potts dynamics on signed CSR graphs -- no file I/O, works with
 * both NX and GT graphs, and is seed-reproducible (unlike the C-subprocess
 * backend, which self-seeds from wall-clock/PID). Reuses the shared
 * potts_metropolis_sweep / calc_potts_energy kernels verbatim, so the Python,
 * C-subprocess and pybind11 backends apply identical update logic (modulo RNG
 * stream: the Python loop draws from numpy, the C kernels from SFMT, so py<->pb
 * agree only in distribution, not bit-for-bit).
 *
 * Naming convention (Python side):
 *   runlang = "pb" / "pb_potts"  -> potts_sampling()
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
#include "LRGSG_potts.h"
#include "sfmtrng.h"
}

namespace py = pybind11;

/* ==================================================================
 * Helper: build NodeEdges adjacency from flat numpy CSR arrays
 * (neigh_indices/weights concatenated, neigh_ptr row pointers). Mirrors
 * the GraphCSR helper in voter_native.cpp.
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
 * Identical scheme to voter_native.cpp::seed_rng. */
static void seed_rng(uint64_t seed_val) {
    uint32_t seed_arr[4];
    seed_arr[0] = static_cast<uint32_t>(seed_val & 0xFFFFFFFFu);
    seed_arr[1] = static_cast<uint32_t>((seed_val >> 32) & 0xFFFFFFFFu);
    seed_arr[2] = static_cast<uint32_t>((seed_val * 0x9E3779B97F4A7C15ULL) & 0xFFFFFFFFu);
    seed_arr[3] = static_cast<uint32_t>((seed_val ^ 0xBE11AC1A0ULL) & 0xFFFFFFFFu) | 1u;
    sfmt_init_by_array(&sfmt, seed_arr, 4);
}

/* Potts order parameter (q * max_frac - 1) / (q - 1); matches
 * PottsModel.order_parameter (q == 1 -> 1.0). */
static double potts_order_parameter(size_t N, const int32_t* s, int q) {
    if (q <= 1) return 1.0;
    std::vector<size_t> counts(static_cast<size_t>(q), 0);
    for (size_t i = 0; i < N; ++i) {
        int32_t v = s[i];
        if (v >= 0 && v < q) counts[static_cast<size_t>(v)]++;
    }
    size_t mx = 0;
    for (size_t c = 0; c < counts.size(); ++c) if (counts[c] > mx) mx = counts[c];
    double max_frac = static_cast<double>(mx) / static_cast<double>(N);
    return (q * max_frac - 1.0) / (q - 1.0);
}

/* ==================================================================
 * Potts sampling: n_sweeps Metropolis sweeps at (q, T).
 *
 * When save_observables is set, the per-state energy (calc_potts_energy) and
 * order parameter are recorded AFTER each sweep -- matching PottsModel.run_py,
 * which appends both inside the sweep loop -- giving n_sweeps records.
 *
 * Returns (final_s[int32], ene[float64], magn[float64]); ene/magn are empty
 * unless save_observables.
 * ================================================================== */
static py::tuple potts_sampling(
    py::array_t<int32_t, py::array::c_style> s_in,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr,
    int q, double T,
    size_t n_sweeps, uint64_t seed_val,
    bool save_observables
) {
    auto s_buf = s_in.request();
    size_t N = static_cast<size_t>(s_buf.size);

    std::vector<int32_t> spins(N);
    std::memcpy(spins.data(), s_buf.ptr, N * sizeof(int32_t));

    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    std::vector<double> ene_out, magn_out;
    if (save_observables) { ene_out.reserve(n_sweeps); magn_out.reserve(n_sweeps); }

    {
        py::gil_scoped_release release;
        seed_rng(seed_val);

        int32_t* s = spins.data();
        size_t* nlen_ptr = graph.nlen.data();
        NodeEdges* ne = graph.node_edges.data();

        for (size_t step = 0; step < n_sweeps; ++step) {
            potts_metropolis_sweep(N, s, q, T, nlen_ptr, ne);
            if (save_observables) {
                ene_out.push_back(calc_potts_energy(N, s, nlen_ptr, ne));
                magn_out.push_back(potts_order_parameter(N, s, q));
            }
        }
    }

    py::array_t<int32_t> out_s(N);
    std::memcpy(out_s.mutable_data(), spins.data(), N * sizeof(int32_t));

    py::array_t<double> out_ene(ene_out.size());
    if (!ene_out.empty())
        std::memcpy(out_ene.mutable_data(), ene_out.data(),
                    ene_out.size() * sizeof(double));
    py::array_t<double> out_magn(magn_out.size());
    if (!magn_out.empty())
        std::memcpy(out_magn.mutable_data(), magn_out.data(),
                    magn_out.size() * sizeof(double));

    return py::make_tuple(out_s, out_ene, out_magn);
}

PYBIND11_MODULE(_potts_native, m) {
    m.doc() = "Native (pybind11) q-state Potts Metropolis kernel on signed CSR graphs.";

    m.def("potts_sampling", &potts_sampling,
        R"pbdoc(
Run q-state Potts Metropolis dynamics on a signed CSR graph.

Parameters
----------
s, neigh_indices, neigh_weights, neigh_ptr : ndarray
    Initial int32 state (values in [0, q)) and the CSR adjacency
    (symmetric, signed weights).
q, T : number of Potts states and temperature.
n_sweeps, seed : run controls (seed-reproducible).
save_observables : record per-state energy + order parameter each sweep.

Returns
-------
tuple[ndarray, ndarray, ndarray]
    (final_state[int32], energy[float64], order_parameter[float64]); the
    two observable arrays have length n_sweeps when save_observables, else 0.
        )pbdoc",
        py::arg("s"), py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("q"), py::arg("T"),
        py::arg("n_sweeps"), py::arg("seed"),
        py::arg("save_observables") = false);
}
