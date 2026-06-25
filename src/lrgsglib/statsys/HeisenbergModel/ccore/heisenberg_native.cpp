/**
 * @file heisenberg_native.cpp
 * @brief Pybind11 wrapper for the 3D Heisenberg Metropolis C kernel.
 *
 * In-process Heisenberg dynamics on signed CSR graphs -- no file I/O, works
 * with both NX and GT graphs, and is seed-reproducible (unlike the C-subprocess
 * backend, which self-seeds from wall-clock/PID). Reuses the shared
 * heisenberg_metropolis_sweep / calc_heisenberg_energy /
 * calc_heisenberg_magnetisation kernels verbatim, so the Python, C-subprocess
 * and pybind11 backends apply identical update logic (modulo RNG stream: the
 * Python loop draws from numpy, the C kernels from SFMT, so py<->pb agree only
 * in distribution, not bit-for-bit).
 *
 * State layout: each node carries a 3D unit vector. Python side holds an
 * (N, 3) float64 array; here it is marshalled as a flat, contiguous 3N double
 * buffer (spin[i*3 + c]) -- exactly the layout the C kernel expects.
 *
 * Naming convention (Python side):
 *   runlang = "pb" / "pb_heisenberg"  -> heisenberg_sampling()
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
#include "LRGSG_heisenberg.h"
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
 * Heisenberg sampling: n_sweeps Metropolis sweeps at (T, delta).
 *
 * State is marshalled as a flat 3N double buffer (spin[i*3 + c]); the input
 * (N, 3) float64 array is required c_style so its memory is already in that
 * layout. When save_observables is set, the per-state energy
 * (calc_heisenberg_energy) and magnetisation (calc_heisenberg_magnetisation)
 * are recorded AFTER each sweep -- matching HeisenbergModel.run_py, which
 * appends both inside the sweep loop -- giving n_sweeps records.
 *
 * Returns (final_s[float64, shape (N,3)], ene[float64], magn[float64]);
 * ene/magn are empty unless save_observables.
 * ================================================================== */
static py::tuple heisenberg_sampling(
    py::array_t<double, py::array::c_style | py::array::forcecast> s_in,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr,
    double T, double delta,
    size_t n_sweeps, uint64_t seed_val,
    bool save_observables
) {
    auto s_buf = s_in.request();
    /* Total scalar count = 3N; recover N from the flattened buffer. */
    size_t total = static_cast<size_t>(s_buf.size);
    size_t N = total / 3;

    std::vector<double> spin(total);
    std::memcpy(spin.data(), s_buf.ptr, total * sizeof(double));

    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    std::vector<double> ene_out, magn_out;
    if (save_observables) { ene_out.reserve(n_sweeps); magn_out.reserve(n_sweeps); }

    {
        py::gil_scoped_release release;
        seed_rng(seed_val);

        double* s = spin.data();
        size_t* nlen_ptr = graph.nlen.data();
        NodeEdges* ne = graph.node_edges.data();

        for (size_t step = 0; step < n_sweeps; ++step) {
            heisenberg_metropolis_sweep(N, s, T, delta, nlen_ptr, ne);
            if (save_observables) {
                ene_out.push_back(calc_heisenberg_energy(N, s, nlen_ptr, ne));
                magn_out.push_back(calc_heisenberg_magnetisation(N, s));
            }
        }
    }

    /* Return final state as (N, 3) float64. */
    py::array_t<double> out_s({static_cast<py::ssize_t>(N),
                               static_cast<py::ssize_t>(3)});
    std::memcpy(out_s.mutable_data(), spin.data(), total * sizeof(double));

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

PYBIND11_MODULE(_heisenberg_native, m) {
    m.doc() = "Native (pybind11) 3D Heisenberg Metropolis kernel on signed CSR graphs.";

    m.def("heisenberg_sampling", &heisenberg_sampling,
        R"pbdoc(
Run 3D Heisenberg Metropolis dynamics on a signed CSR graph.

Parameters
----------
s : ndarray, shape (N, 3), float64
    Initial state of unit vectors (one 3D spin per node).
neigh_indices, neigh_weights, neigh_ptr : ndarray
    The CSR adjacency (symmetric, signed weights).
T, delta : temperature and maximum proposal rotation angle.
n_sweeps, seed : run controls (seed-reproducible).
save_observables : record per-state energy + magnetisation each sweep.

Returns
-------
tuple[ndarray, ndarray, ndarray]
    (final_state[float64, (N,3)], energy[float64], magnetisation[float64]);
    the two observable arrays have length n_sweeps when save_observables,
    else 0.
        )pbdoc",
        py::arg("s"), py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("T"), py::arg("delta"),
        py::arg("n_sweeps"), py::arg("seed"),
        py::arg("save_observables") = false);
}
