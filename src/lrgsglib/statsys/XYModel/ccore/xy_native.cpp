/**
 * @file xy_native.cpp
 * @brief Pybind11 wrapper for the XY (planar rotator) Metropolis C kernel.
 *
 * In-process XY dynamics on signed CSR graphs -- no file I/O, works with both
 * NX and GT graphs, and is seed-reproducible (unlike the C-subprocess backend,
 * which self-seeds from wall-clock/PID). Reuses the shared xy_metropolis_sweep
 * / calc_xy_energy / calc_xy_magnetisation kernels verbatim, so the Python,
 * C-subprocess and pybind11 backends apply identical update logic (modulo RNG
 * stream: the Python loop draws from numpy, the C kernels from SFMT, so py<->pb
 * agree only in distribution, not bit-for-bit).
 *
 * Naming convention (Python side):
 *   runlang = "pb" / "pb_xy"  -> xy_sampling()
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
#include "LRGSG_xy.h"
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
 * XY sampling: n_sweeps Metropolis sweeps at (T, delta).
 *
 * When save_observables is set, the per-state energy (calc_xy_energy) and
 * magnetisation (calc_xy_magnetisation) are recorded AFTER each sweep --
 * matching XYModel.run_py, which appends both inside the sweep loop -- giving
 * n_sweeps records.
 *
 * Returns (final_theta[float64], ene[float64], magn[float64]); ene/magn are
 * empty unless save_observables.
 * ================================================================== */
static py::tuple xy_sampling(
    py::array_t<double, py::array::c_style> theta_in,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr,
    double T, double delta,
    size_t n_sweeps, uint64_t seed_val,
    bool save_observables
) {
    auto t_buf = theta_in.request();
    size_t N = static_cast<size_t>(t_buf.size);

    std::vector<double> theta(N);
    std::memcpy(theta.data(), t_buf.ptr, N * sizeof(double));

    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    std::vector<double> ene_out, magn_out;
    if (save_observables) { ene_out.reserve(n_sweeps); magn_out.reserve(n_sweeps); }

    {
        py::gil_scoped_release release;
        seed_rng(seed_val);

        double* th = theta.data();
        size_t* nlen_ptr = graph.nlen.data();
        NodeEdges* ne = graph.node_edges.data();

        for (size_t step = 0; step < n_sweeps; ++step) {
            xy_metropolis_sweep(N, th, T, delta, nlen_ptr, ne);
            if (save_observables) {
                ene_out.push_back(calc_xy_energy(N, th, nlen_ptr, ne));
                magn_out.push_back(calc_xy_magnetisation(N, th));
            }
        }
    }

    py::array_t<double> out_theta(N);
    std::memcpy(out_theta.mutable_data(), theta.data(), N * sizeof(double));

    py::array_t<double> out_ene(ene_out.size());
    if (!ene_out.empty())
        std::memcpy(out_ene.mutable_data(), ene_out.data(),
                    ene_out.size() * sizeof(double));
    py::array_t<double> out_magn(magn_out.size());
    if (!magn_out.empty())
        std::memcpy(out_magn.mutable_data(), magn_out.data(),
                    magn_out.size() * sizeof(double));

    return py::make_tuple(out_theta, out_ene, out_magn);
}

PYBIND11_MODULE(_xy_native, m) {
    m.doc() = "Native (pybind11) XY (planar rotator) Metropolis kernel on signed CSR graphs.";

    m.def("xy_sampling", &xy_sampling,
        R"pbdoc(
Run XY (planar rotator) Metropolis dynamics on a signed CSR graph.

Parameters
----------
theta, neigh_indices, neigh_weights, neigh_ptr : ndarray
    Initial float64 angle state (radians) and the CSR adjacency
    (symmetric, signed weights).
T, delta : temperature and maximum angular perturbation (radians).
n_sweeps, seed : run controls (seed-reproducible).
save_observables : record per-state energy + magnetisation each sweep.

Returns
-------
tuple[ndarray, ndarray, ndarray]
    (final_theta[float64], energy[float64], magnetisation[float64]); the
    two observable arrays have length n_sweeps when save_observables, else 0.
        )pbdoc",
        py::arg("theta"), py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("T"), py::arg("delta"),
        py::arg("n_sweeps"), py::arg("seed"),
        py::arg("save_observables") = false);
}
