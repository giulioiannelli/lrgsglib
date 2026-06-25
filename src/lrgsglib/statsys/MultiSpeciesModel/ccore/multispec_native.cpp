/**
 * @file multispec_native.cpp
 * @brief Pybind11 wrapper for the multi-species Metropolis C kernel.
 *
 * In-process multi-species (k components/node, q states) dynamics on signed
 * CSR graphs -- no file I/O, works with both NX and GT graphs, and is
 * seed-reproducible (unlike the C-subprocess backend, which self-seeds from
 * wall-clock/PID). Reuses the shared ``multispec_metropolis_sweep`` kernel
 * verbatim, so the C-subprocess and pybind11 backends apply identical update
 * logic (identity interaction matrix, uniform q across species).
 *
 * State layout: row-major int32 of length N*k, s[nd*k + comp] -- exactly the
 * flattened (N, k) array marshalled by MultiSpeciesModel.
 *
 * Naming convention (Python side):
 *   runlang = "pb" / "pb_multispec"  -> multispec_sampling()
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
#include "LRGSG_vecdynsys.h"
#include "LRGSG_multispec.h"
#include "sfmtrng.h"
}

namespace py = pybind11;

/* ==================================================================
 * Helper: build NodeEdges adjacency from flat numpy CSR arrays
 * (neigh_indices/weights concatenated, neigh_ptr row pointers). Mirrors
 * the GraphCSR helper in potts_native.cpp.
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

/* Multi-species energy with identity interaction matrix + uniform q:
 *   H = -sum_{(i,j)} sum_s A_ij * delta(s_i^s, s_j^s)
 * Matches MultiSpeciesModel.energy() in the C-backend regime (identity
 * interaction, only same-component coupling contributes). Each undirected
 * edge appears twice in the CSR adjacency, so the half-sum below recovers the
 * (i<j) Python convention. */
static double calc_multispec_energy(size_t N, const int32_t* s, int k,
                                    const size_t* nlen,
                                    const NodeEdges* ne) {
    double E = 0.0;
    for (size_t i = 0; i < N; ++i) {
        size_t deg = nlen[i];
        for (size_t e = 0; e < deg; ++e) {
            size_t j = ne[i].neighbors[e];
            double w = ne[i].weights[e];
            for (int comp = 0; comp < k; ++comp) {
                E -= w * kronecker_delta_i32(s[i * k + comp], s[j * k + comp]);
            }
        }
    }
    return 0.5 * E;
}

/* ==================================================================
 * Multi-species sampling: n_sweeps Metropolis sweeps at (k, q, T).
 *
 * When save_observables is set, the per-state energy is recorded AFTER each
 * sweep -- matching MultiSpeciesModel.run_py, which appends inside the sweep
 * loop -- giving n_sweeps records.
 *
 * Returns (final_s[int32, length N*k], ene[float64]); ene is empty unless
 * save_observables.
 * ================================================================== */
static py::tuple multispec_sampling(
    py::array_t<int32_t, py::array::c_style> s_in,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr,
    size_t N, int k, int q, double T,
    size_t n_sweeps, uint64_t seed_val,
    bool save_observables
) {
    auto s_buf = s_in.request();
    size_t total = static_cast<size_t>(s_buf.size);  // == N * k

    std::vector<int32_t> spins(total);
    std::memcpy(spins.data(), s_buf.ptr, total * sizeof(int32_t));

    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    std::vector<double> ene_out;
    if (save_observables) ene_out.reserve(n_sweeps);

    {
        py::gil_scoped_release release;
        seed_rng(seed_val);

        int32_t* s = spins.data();
        size_t* nlen_ptr = graph.nlen.data();
        NodeEdges* ne = graph.node_edges.data();

        for (size_t step = 0; step < n_sweeps; ++step) {
            multispec_metropolis_sweep(N, s, k, q, T, nlen_ptr, ne);
            if (save_observables) {
                ene_out.push_back(calc_multispec_energy(N, s, k, nlen_ptr, ne));
            }
        }
    }

    py::array_t<int32_t> out_s(total);
    std::memcpy(out_s.mutable_data(), spins.data(), total * sizeof(int32_t));

    py::array_t<double> out_ene(ene_out.size());
    if (!ene_out.empty())
        std::memcpy(out_ene.mutable_data(), ene_out.data(),
                    ene_out.size() * sizeof(double));

    return py::make_tuple(out_s, out_ene);
}

PYBIND11_MODULE(_multispec_native, m) {
    m.doc() = "Native (pybind11) multi-species Metropolis kernel on signed CSR graphs.";

    m.def("multispec_sampling", &multispec_sampling,
        R"pbdoc(
Run multi-species Metropolis dynamics on a signed CSR graph.

Parameters
----------
s, neigh_indices, neigh_weights, neigh_ptr : ndarray
    Initial int32 state (flattened (N, k), values in [0, q)) and the CSR
    adjacency (symmetric, signed weights).
N, k, q, T : node count, species count, states per species, temperature.
n_sweeps, seed : run controls (seed-reproducible).
save_observables : record per-state energy each sweep.

Returns
-------
tuple[ndarray, ndarray]
    (final_state[int32, length N*k], energy[float64]); the energy array has
    length n_sweeps when save_observables, else 0.
        )pbdoc",
        py::arg("s"), py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("N"), py::arg("k"), py::arg("q"),
        py::arg("T"), py::arg("n_sweeps"), py::arg("seed"),
        py::arg("save_observables") = false);
}
