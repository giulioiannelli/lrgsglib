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
#include <cstdlib>
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
    bool savemagn,
    int rule, size_t q, double eps, double alpha,
    int upd_mode, bool absorbing,
    bool track_clusters, int cluster_mode, bool savedyn,
    const py::array_t<int64_t>& snap_sweeps
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
    if (savemagn) magn_out.reserve(n_sweeps);
    long absorbed_at = -1;

    std::vector<int8_t> snew(mode == VOTER_UPD_SYNC ? N : 0);
    std::vector<size_t> cdeg;
    size_t total = 0;

    /* Cluster-size-distribution capture (flattened; built into py objects below,
     * with the GIL held). Only populated when track_clusters is requested. */
    std::vector<size_t> cl_off, cl_size, cl_cnt;
    size_t cl_nrec = 0;

    /* Optional per-sweep trajectory (savedyn): row-major (n_rec, N) int8. */
    std::vector<int8_t> snap_buf;
    size_t snap_nrec = 0;
    /* Optional sorted sweep indices to record (log-spaced sampling); empty = all. */
    std::vector<size_t> snap_sweeps_vec;
    {
        auto ss = snap_sweeps.unchecked<1>();
        snap_sweeps_vec.reserve(static_cast<size_t>(ss.size()));
        for (py::ssize_t k = 0; k < ss.size(); ++k)
            snap_sweeps_vec.push_back(static_cast<size_t>(ss(k)));
    }

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

        /* optional cluster recompute, shared by every schedule (NULL => skipped).
         * Created with s = spins.data(); the non-gillespie loop re-points it at the
         * live buffer each record via clusters_set_state (sync swaps s). */
        ClusterCtx *cctx = track_clusters
            ? clusters_create(N, s, nlen_ptr, ne, cluster_mode == 1 ? 1 : 0)
            : nullptr;

        if (mode == VOTER_UPD_GILLESPIE) {
            /* Rejection-free CTMC (shared _ccore kernel); magn sampled at integer
             * sweep times, t_run shortens if it freezes (absorbing). */
            std::vector<double> magn_buf(savemagn ? n_sweeps : 0);
            double *magn_ptr = savemagn ? magn_buf.data() : nullptr;
            /* optional full trajectory (savedyn): the kernel grows a buffer to the
             * ACTUAL recorded-sweep count (never the full n_sweeps*N) and hands it
             * back here for us to copy out and free. */
            spin_tp snap_raw = nullptr;
            spin_tp *snap_out = savedyn ? &snap_raw : nullptr;
            const size_t *ss_ptr =
                snap_sweeps_vec.empty() ? nullptr : snap_sweeps_vec.data();
            size_t snap_rows = 0;
            size_t t_run = voter_ctmc_run(
                N, s, nlen_ptr, ne, n_sweeps,
                savemagn ? 1 : 0, magn_ptr, snap_out, ss_ptr,
                snap_sweeps_vec.size(), &snap_rows,
                absorbing ? 1 : 0, &absorbed_at, cctx);
            if (savemagn)
                magn_out.assign(magn_buf.begin(), magn_buf.begin() + t_run);
            if (savedyn && snap_raw) {
                snap_nrec = snap_rows;
                snap_buf.assign(snap_raw, snap_raw + snap_rows * N);
                free(snap_raw);
            }
        } else {
            size_t snap_si = 0;
            for (size_t step = 0; step < n_sweeps; ++step) {
                if (savemagn) magn_out.push_back(calc_magn(N, s));
                if (savedyn) {
                    bool take = snap_sweeps_vec.empty() ||
                        (snap_si < snap_sweeps_vec.size() &&
                         snap_sweeps_vec[snap_si] == step);
                    if (take) {
                        snap_buf.insert(snap_buf.end(), s, s + N);
                        ++snap_nrec;
                        if (!snap_sweeps_vec.empty()) ++snap_si;
                    }
                }
                if (cctx) { clusters_set_state(cctx, s); clusters_record(cctx); }
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
        if (cctx) {                       /* copy out raw cluster records, then free */
            cl_nrec = clusters_num_records(cctx);
            const size_t *off = clusters_offsets(cctx);
            const size_t *sz  = clusters_sizes(cctx);
            const size_t *ct  = clusters_counts(cctx);
            size_t nent = off[cl_nrec];
            cl_off.assign(off, off + cl_nrec + 1);
            cl_size.assign(sz, sz + nent);
            cl_cnt.assign(ct, ct + nent);
            clusters_free(cctx);
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

    /* cluster-size distribution: one dict {size: count} per recorded sweep
     * (empty list unless track_clusters was requested). Built with the GIL held. */
    py::list cluster_dist;
    for (size_t r = 0; r < cl_nrec; ++r) {
        py::dict d;
        for (size_t k = cl_off[r]; k < cl_off[r + 1]; ++k)
            d[py::int_(cl_size[k])] = py::int_(cl_cnt[k]);
        cluster_dist.append(d);
    }

    /* full per-sweep trajectory (savedyn): (n_rec, N) int8 array, empty if off */
    std::vector<py::ssize_t> s_t_shape{(py::ssize_t)snap_nrec, (py::ssize_t)N};
    py::array_t<int8_t> out_s_t(s_t_shape);
    if (snap_nrec)
        std::memcpy(out_s_t.mutable_data(), snap_buf.data(),
                    (size_t)snap_nrec * N * sizeof(int8_t));

    return py::make_tuple(out_spins, out_magn, absorbed_at, cluster_dist, out_s_t);
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
n_sweeps, seed, savemagn : run controls.
rule : 0 linear | 1 majority | 2 qvoter | 3 nonlinear.
q, eps, alpha : q-voter / nonlinear parameters.
upd_mode : 0 async | 1 sync | 2 link | 3 gillespie (rejection-free CTMC, linear).
absorbing : stop at the first zero-frustration configuration.
track_clusters : record the cluster-size distribution per sweep (any schedule;
    full connected-components recompute over active edges at each record).
cluster_mode : 0 satisfied (signed domains) | 1 rawspin (edge sign ignored).
savedyn : record the per-sweep spin trajectory (returned as the 5th element).
snap_sweeps : optional sorted sweep indices to record (log-spaced sampling);
    empty array = every sweep. Bounds the in-RAM trajectory to len(snap_sweeps) rows.

Returns
-------
tuple[ndarray, ndarray, int, list, ndarray]
    (final_spins[int8], magnetization[float64] of length = sweeps run,
     absorbed_at, cluster_dist, s_t) where absorbed_at = -1 if it never froze,
     cluster_dist is a list (one {size: count} dict per recorded sweep; empty
     unless track_clusters), and s_t is an (n_recorded, N) int8 array of the
     per-sweep configurations (shape (0, N) unless savedyn).
        )pbdoc",
        py::arg("spins"), py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("n_sweeps"), py::arg("seed"),
        py::arg("savemagn") = true,
        py::arg("rule") = 0, py::arg("q") = 2, py::arg("eps") = 0.0,
        py::arg("alpha") = 1.0, py::arg("upd_mode") = 0,
        py::arg("absorbing") = false,
        py::arg("track_clusters") = false, py::arg("cluster_mode") = 0,
        py::arg("savedyn") = false,
        py::arg("snap_sweeps") = py::array_t<int64_t>());
}
