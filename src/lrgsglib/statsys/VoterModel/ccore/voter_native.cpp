/**
 * @file voter_native.cpp
 * @brief Pybind11 wrapper for the voter-model C kernels.
 *
 * In-process voter dynamics on signed CSR graphs -- no file I/O, works with
 * both NX and GT graphs, and is seed-reproducible (unlike the C-subprocess
 * backend, which self-seeds from wall-clock/PID). Reuses the shared
 * voter_apply_rule / voter_*_step kernels, so the Python, C-subprocess and
 * pybind11 backends apply identical update logic across the Phase-3 rule family
 * (Axis A) and sampler axis (Axis B), with optional absorbing-state early stop.
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
#include "LRGSG_ctmc.h"
#include "LRGSG_utils.h"
#include "LRGSG_vm.h"
#include "sfmtrng.h"
}

namespace py = pybind11;

/* ==================================================================
 * Helper: build NodesEdges adjacency from flat numpy CSR arrays.
 *
 *   neigh_indices : int64[total_degree]  (concatenated neighbour lists)
 *   neigh_weights : float64[total_degree](concatenated signed weights)
 *   neigh_ptr     : int64[N+1]           (row pointers into the above)
 * ================================================================== */
struct GraphCSR {
    size_t N;
    std::vector<size_t> nlen;           /* degree of each node */
    std::vector<NodeEdges> node_edges;  /* per-node adjacency  */

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

/* Seed the global SFMT state directly from seed_val (reproducible in-process
 * runs; __set_seed_SFMT() is deliberately not used -- it self-seeds). */
static void seed_rng(uint64_t seed_val) {
    uint32_t seed_arr[4];
    seed_arr[0] = static_cast<uint32_t>(seed_val & 0xFFFFFFFFu);
    seed_arr[1] = static_cast<uint32_t>((seed_val >> 32) & 0xFFFFFFFFu);
    seed_arr[2] = static_cast<uint32_t>((seed_val * 0x9E3779B97F4A7C15ULL) & 0xFFFFFFFFu);
    seed_arr[3] = static_cast<uint32_t>((seed_val ^ 0xBE11AC1A0ULL) & 0xFFFFFFFFu) | 1u;
    sfmt_init_by_array(&sfmt, seed_arr, 4);
}

/* ==================================================================
 * Voter sampling: rule family (Axis A) x sampler (Axis B) + absorbing stop.
 *
 * Runs up to `n_sweeps` sweeps under `rule`/`upd_mode`, recording the
 * magnetization before each sweep (when requested). With `absorbing=true` the
 * run stops at the first zero-frustration configuration.
 *
 * Returns (final_spins[int8], magn[float64] of length = sweeps run,
 *          absorbed_at[int]) where absorbed_at = -1 if it never froze.
 * ================================================================== */
static py::tuple voter_sampling(
    py::array_t<int8_t, py::array::c_style> spins_in,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr,
    size_t n_sweeps,
    uint64_t seed_val,
    bool save_magnetization,
    int rule, size_t q, double eps, double alpha,
    int upd_mode, bool absorbing
) {
    auto s_buf = spins_in.request();
    size_t N = static_cast<size_t>(s_buf.size);

    std::vector<int8_t> spins(N);
    std::memcpy(spins.data(), s_buf.ptr, N * sizeof(int8_t));

    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    voter_params vp;
    vp.rule  = static_cast<voter_rule_t>(rule);
    vp.q     = q;
    vp.eps   = eps;
    vp.alpha = alpha;
    voter_upd_t mode = static_cast<voter_upd_t>(upd_mode);

    std::vector<double> magn_out;
    if (save_magnetization) magn_out.reserve(n_sweeps);
    long absorbed_at = -1;

    std::vector<int8_t> snew(mode == VOTER_UPD_SYNC ? N : 0);
    std::vector<size_t> cdeg;
    size_t total = 0;

    {
        py::gil_scoped_release release;
        seed_rng(seed_val);

        spin_tp s = spins.data();
        size_tp nlen_ptr = graph.nlen.data();
        NodesEdges ne = graph.node_edges.data();

        if (mode == VOTER_UPD_LINK) {
            cdeg.resize(N + 1);
            total = voter_build_cdeg(N, nlen_ptr, cdeg.data());
        }
        spin_tp sbuf = snew.empty() ? nullptr : snew.data();

        if (mode == VOTER_UPD_GILLESPIE) {
            /* Rejection-free CTMC (shared _ccore kernel); magn sampled at integer
             * sweep times, t_run shortens if it freezes (absorbing). */
            std::vector<double> magn_buf(save_magnetization ? n_sweeps : 0);
            double *magn_ptr = save_magnetization ? magn_buf.data() : nullptr;
            size_t t_run = voter_ctmc_run(
                N, s, nlen_ptr, ne, n_sweeps,
                save_magnetization ? 1 : 0, magn_ptr,
                absorbing ? 1 : 0, &absorbed_at);
            if (save_magnetization)
                magn_out.assign(magn_buf.begin(), magn_buf.begin() + t_run);
        } else {
            for (size_t step = 0; step < n_sweeps; ++step) {
                if (save_magnetization) magn_out.push_back(calc_magn(N, s));
                if (absorbing && voter_count_frustrated(N, s, nlen_ptr, ne) == 0) {
                    absorbed_at = static_cast<long>(step);
                    break;
                }
                if (mode == VOTER_UPD_SYNC) {
                    voter_sync_step(N, s, sbuf, nlen_ptr, ne, vp);
                    spin_tp tmp = s; s = sbuf; sbuf = tmp;   /* swap */
                } else if (mode == VOTER_UPD_LINK) {
                    voter_link_step(N, s, nlen_ptr, ne, cdeg.data(), total);
                } else {
                    voter_model_Nstep(N, s, nlen_ptr, ne, vp);
                }
            }
        }
        /* `s` may point at `snew`'s storage after an odd number of sync swaps;
         * copy the live state back into `spins` for the return value. */
        if (s != spins.data())
            std::memcpy(spins.data(), s, N * sizeof(int8_t));
    }

    py::array_t<int8_t> out_spins(N);
    std::memcpy(out_spins.mutable_data(), spins.data(), N * sizeof(int8_t));

    py::array_t<double> out_magn(magn_out.size());
    if (!magn_out.empty())
        std::memcpy(out_magn.mutable_data(), magn_out.data(),
                    magn_out.size() * sizeof(double));

    return py::make_tuple(out_spins, out_magn, absorbed_at);
}

PYBIND11_MODULE(_voter_native, m) {
    m.doc() = "Native (pybind11) voter-model kernels on signed CSR graphs.";

    m.def("voter_sampling", &voter_sampling,
        R"pbdoc(
Run voter dynamics on a signed CSR graph (rule family x sampler axis).

Parameters
----------
spins, neigh_indices, neigh_weights, neigh_ptr : ndarray
    Initial int8 state and the CSR adjacency (symmetric, signed weights).
n_sweeps, seed, save_magnetization : run controls.
rule : 0 linear | 1 majority | 2 qvoter | 3 nonlinear.
q, eps, alpha : q-voter / nonlinear parameters.
upd_mode : 0 async | 1 sync | 2 link | 3 gillespie (rejection-free CTMC, linear).
absorbing : stop at the first zero-frustration configuration.

Returns
-------
tuple[ndarray, ndarray, int]
    (final_spins[int8], magnetization[float64] of length = sweeps run,
     absorbed_at) where absorbed_at = -1 if it never froze.
        )pbdoc",
        py::arg("spins"), py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("n_sweeps"), py::arg("seed"),
        py::arg("save_magnetization") = true,
        py::arg("rule") = 0, py::arg("q") = 2, py::arg("eps") = 0.0,
        py::arg("alpha") = 1.0, py::arg("upd_mode") = 0,
        py::arg("absorbing") = false);
}
