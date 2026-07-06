/**
 * @file ising_native.cpp
 * @brief Pybind11 wrapper for Ising dynamics C kernels.
 *
 * Provides Python-callable functions that wrap the existing C Metropolis,
 * Simulated Annealing, and Parallel Tempering kernels.  Data is passed
 * directly via numpy arrays -- no file I/O is needed.
 *
 * Naming convention (Python side):
 *   runlang = "pb_met"  → metropolis_sampling()
 *   runlang = "pb_sa"   → sa_sampling()
 *   runlang = "pb_pt"   → pt_sampling()
 */
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <vector>

/* ------------------------------------------------------------------ */
/* Pull in the existing C kernels                                      */
/* ------------------------------------------------------------------ */
extern "C" {
#include "LRGSG_rbim.h"
#include "LRGSG_sa.h"
#include "LRGSG_pt.h"
#include "LRGSG_cluster.h"
#include "LRGSG_utils.h"
#include "sfmtrng.h"
}

namespace py = pybind11;

/* ==================================================================
 * Helper: build NodesEdges from flat numpy arrays
 * ==================================================================
 *
 * The C kernels expect a `NodesEdges` structure (one NodeEdges per
 * node, each containing a `neighbors` array and a `weights` array).
 *
 * Python passes the graph as CSR-like arrays:
 *   neigh_indices : int64[total_degree]    (concatenated neighbor lists)
 *   neigh_weights : float64[total_degree]  (concatenated weight lists)
 *   neigh_ptr     : int64[N+1]             (row pointers into above)
 *   neigh_len     : int64[N]               (degree of each node)
 */

struct GraphCSR {
    size_t N;
    std::vector<size_t> nlen;           /* degree of each node */
    std::vector<NodeEdges> node_edges;  /* per-node adjacency */

    /* Internal storage that the NodeEdges point into */
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

/* ==================================================================
 * Helper: seed the SFMT RNG
 * ================================================================== */
static void seed_rng(uint64_t seed_val) {
    static uint32_t seed_arr[4];
    seed_arr[0] = static_cast<uint32_t>(seed_val & 0xFFFFFFFF);
    seed_arr[1] = static_cast<uint32_t>((seed_val >> 32) & 0xFFFFFFFF);
    seed_arr[2] = static_cast<uint32_t>(seed_val * 0x9E3779B97F4A7C15ULL);
    seed_arr[3] = static_cast<uint32_t>(seed_val ^ 0xBE11AC1A0ULL);
    seed_rand = seed_arr;
    __set_seed_SFMT();
    /* The asynchronous visit order (__gen_rand_u64_array in sfmtrng.c)
     * derives each per-sweep sub-seed from libc rand(); seed that stream
     * too, or seeded runs do not reproduce their site sequences. */
    srand(static_cast<unsigned>(seed_val & 0xFFFFFFFF));
}

/* ==================================================================
 * Metropolis sampling
 * ==================================================================
 *
 * Runs `n_sweeps` MC sweeps at temperature `T`, recording energy and
 * magnetization after each sweep.  Returns (final_spins, energy, magn).
 */
static py::tuple metropolis_sampling(
    py::array_t<int8_t, py::array::c_style> spins_in,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr,
    const py::array_t<double>&  field,
    double T,
    size_t n_sweeps,
    uint64_t seed_val,
    const std::string& update_mode
) {
    auto s_buf = spins_in.request();
    size_t N = static_cast<size_t>(s_buf.size);

    /* Copy input spins so the original is not mutated */
    std::vector<int8_t> spins(N);
    std::memcpy(spins.data(), s_buf.ptr, N * sizeof(int8_t));

    /* Build graph CSR */
    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    /* Get field pointer */
    auto h_buf = field.request();
    const double *h_ptr = static_cast<const double*>(h_buf.ptr);

    /* Output arrays */
    std::vector<double> ene_out(n_sweeps);
    std::vector<double> magn_out(n_sweeps);

    {
        py::gil_scoped_release release;

        seed_rng(seed_val);
        initialize_glauberMetropolis(T, h_ptr);

        spin_tp s = spins.data();
        size_tp nlen_ptr = graph.nlen.data();
        NodesEdges ne = graph.node_edges.data();
        const char *mode_c = update_mode.c_str();

        for (size_t step = 0; step < n_sweeps; ++step) {
            /* Record observables before sweep */
            ene_out[step]  = calc_totEnergy(N, s, nlen_ptr, ne);
            magn_out[step] = calc_magn(N, s);

            /* MC sweep */
            glauberMetropolis_Nstep(N, T, s, nlen_ptr, ne, mode_c);
        }
    }

    /* Build return arrays */
    py::array_t<int8_t> out_spins(N);
    std::memcpy(out_spins.mutable_data(), spins.data(), N * sizeof(int8_t));

    py::array_t<double> out_ene(n_sweeps);
    std::memcpy(out_ene.mutable_data(), ene_out.data(), n_sweeps * sizeof(double));

    py::array_t<double> out_magn(n_sweeps);
    std::memcpy(out_magn.mutable_data(), magn_out.data(), n_sweeps * sizeof(double));

    return py::make_tuple(out_spins, out_ene, out_magn);
}

/* ==================================================================
 * Simulated Annealing
 * ================================================================== */
static py::tuple sa_sampling(
    py::array_t<int8_t, py::array::c_style> spins_in,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr,
    const py::array_t<double>&  field,
    double T_init,
    double T_final,
    const std::string& schedule_str,
    double cooling_rate,
    size_t steps_per_T,
    size_t n_temps,
    uint64_t seed_val,
    const std::string& update_mode
) {
    auto s_buf = spins_in.request();
    size_t N = static_cast<size_t>(s_buf.size);

    std::vector<int8_t> spins(N);
    std::memcpy(spins.data(), s_buf.ptr, N * sizeof(int8_t));

    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    auto h_buf = field.request();
    const double *h_ptr = static_cast<const double*>(h_buf.ptr);

    size_t total_steps = n_temps * steps_per_T;
    std::vector<double> ene_out(total_steps);
    std::vector<double> magn_out(total_steps);
    std::vector<double> T_out(n_temps);

    sa_cooling_t schedule = sa_parse_schedule(schedule_str.c_str());

    size_t actual_steps;
    {
        py::gil_scoped_release release;

        seed_rng(seed_val);

        spin_tp s = spins.data();
        size_tp nlen_ptr = graph.nlen.data();
        NodesEdges ne = graph.node_edges.data();
        const char *mode_c = update_mode.c_str();

        actual_steps = sa_run(
            N, s, nlen_ptr, ne,
            T_init, T_final, schedule, cooling_rate,
            steps_per_T, n_temps,
            h_ptr, mode_c,
            ene_out.data(), magn_out.data(), T_out.data(),
            0, nullptr  /* no snapshots */
        );
    }

    /* Trim to actual length */
    py::array_t<int8_t> out_spins(N);
    std::memcpy(out_spins.mutable_data(), spins.data(), N * sizeof(int8_t));

    py::array_t<double> out_ene(actual_steps);
    std::memcpy(out_ene.mutable_data(), ene_out.data(), actual_steps * sizeof(double));

    py::array_t<double> out_magn(actual_steps);
    std::memcpy(out_magn.mutable_data(), magn_out.data(), actual_steps * sizeof(double));

    py::array_t<double> out_T(n_temps);
    std::memcpy(out_T.mutable_data(), T_out.data(), n_temps * sizeof(double));

    return py::make_tuple(out_spins, out_ene, out_magn, out_T);
}

/* ==================================================================
 * Parallel Tempering
 * ================================================================== */
static py::tuple pt_sampling(
    py::array_t<int8_t, py::array::c_style> spins_in,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr,
    const py::array_t<double>&  field,
    const py::array_t<double>&  T_ladder_in,
    size_t steps_per_exchange,
    size_t n_exchanges,
    uint64_t seed_val,
    const std::string& update_mode
) {
    auto s_buf = spins_in.request();
    size_t N = static_cast<size_t>(s_buf.size);

    auto tl_buf = T_ladder_in.request();
    size_t n_replicas = static_cast<size_t>(tl_buf.size);

    /* Copy temperature ladder */
    std::vector<double> T_ladder(n_replicas);
    std::memcpy(T_ladder.data(), tl_buf.ptr, n_replicas * sizeof(double));

    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    auto h_buf = field.request();
    const double *h_ptr = static_cast<const double*>(h_buf.ptr);

    /* Initialise replicas (each gets a copy of input spins) */
    std::vector<std::vector<int8_t>> replica_storage(n_replicas, std::vector<int8_t>(N));
    std::vector<spin_tp> replicas(n_replicas);
    for (size_t r = 0; r < n_replicas; ++r) {
        std::memcpy(replica_storage[r].data(), s_buf.ptr, N * sizeof(int8_t));
        replicas[r] = replica_storage[r].data();
    }

    /* Output arrays */
    std::vector<double> ene_out(n_replicas * n_exchanges);
    std::vector<double> magn_out(n_replicas * n_exchanges);
    std::vector<int> exchange_out(n_exchanges * (n_replicas - 1), 0);

    {
        py::gil_scoped_release release;

        seed_rng(seed_val);

        size_tp nlen_ptr = graph.nlen.data();
        NodesEdges ne = graph.node_edges.data();
        const char *mode_c = update_mode.c_str();

        pt_run(
            N, n_replicas, replicas.data(), T_ladder.data(),
            steps_per_exchange, n_exchanges,
            nlen_ptr, ne,
            h_ptr, mode_c,
            ene_out.data(), magn_out.data(), exchange_out.data(),
            0, nullptr  /* no snapshots */
        );
    }

    /* Build output: spins of coldest replica, energy[n_replicas, n_exchanges],
     * magnetization[n_replicas, n_exchanges], exchange_history */
    py::array_t<int8_t> out_spins(N);
    std::memcpy(out_spins.mutable_data(), replicas[0], N * sizeof(int8_t));

    py::array_t<double> out_ene({n_replicas, n_exchanges});
    std::memcpy(out_ene.mutable_data(), ene_out.data(),
                n_replicas * n_exchanges * sizeof(double));

    py::array_t<double> out_magn({n_replicas, n_exchanges});
    std::memcpy(out_magn.mutable_data(), magn_out.data(),
                n_replicas * n_exchanges * sizeof(double));

    py::array_t<int> out_exch({n_exchanges, n_replicas - 1});
    std::memcpy(out_exch.mutable_data(), exchange_out.data(),
                n_exchanges * (n_replicas - 1) * sizeof(int));

    return py::make_tuple(out_spins, out_ene, out_magn, out_exch);
}

/* ==================================================================
 * Single energy calculation (utility)
 * ================================================================== */
/* ==================================================================
 * Wolff cluster sampling
 * ==================================================================
 *
 * Runs `n_sweeps` Wolff sweeps at temperature `T`, recording energy,
 * magnetization, and cluster sizes.
 */
static py::tuple wolff_sampling(
    py::array_t<int8_t, py::array::c_style> spins_in,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr,
    const py::array_t<double>&  field,
    double T,
    size_t n_sweeps,
    uint64_t seed_val
) {
    auto s_buf = spins_in.request();
    size_t N = static_cast<size_t>(s_buf.size);

    std::vector<int8_t> spins(N);
    std::memcpy(spins.data(), s_buf.ptr, N * sizeof(int8_t));

    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    std::vector<double> ene_out(n_sweeps);
    std::vector<double> magn_out(n_sweeps);
    std::vector<int64_t> cluster_out(n_sweeps);

    {
        py::gil_scoped_release release;

        seed_rng(seed_val);

        spin_tp s = spins.data();
        size_tp nlen_ptr = graph.nlen.data();
        NodesEdges ne = graph.node_edges.data();

        for (size_t step = 0; step < n_sweeps; ++step) {
            ene_out[step]  = calc_totEnergy(N, s, nlen_ptr, ne);
            magn_out[step] = calc_magn(N, s);
            cluster_out[step] = static_cast<int64_t>(
                wolff_sweep(N, T, s, nlen_ptr, ne));
        }
    }

    py::array_t<int8_t> out_spins(N);
    std::memcpy(out_spins.mutable_data(), spins.data(), N * sizeof(int8_t));

    py::array_t<double> out_ene(n_sweeps);
    std::memcpy(out_ene.mutable_data(), ene_out.data(), n_sweeps * sizeof(double));

    py::array_t<double> out_magn(n_sweeps);
    std::memcpy(out_magn.mutable_data(), magn_out.data(), n_sweeps * sizeof(double));

    py::array_t<int64_t> out_clusters(n_sweeps);
    std::memcpy(out_clusters.mutable_data(), cluster_out.data(), n_sweeps * sizeof(int64_t));

    return py::make_tuple(out_spins, out_ene, out_magn, out_clusters);
}


/* ==================================================================
 * Swendsen-Wang cluster sampling
 * ================================================================== */
static py::tuple sw_sampling(
    py::array_t<int8_t, py::array::c_style> spins_in,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr,
    const py::array_t<double>&  field,
    double T,
    size_t n_sweeps,
    uint64_t seed_val
) {
    auto s_buf = spins_in.request();
    size_t N = static_cast<size_t>(s_buf.size);

    std::vector<int8_t> spins(N);
    std::memcpy(spins.data(), s_buf.ptr, N * sizeof(int8_t));

    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    std::vector<double> ene_out(n_sweeps);
    std::vector<double> magn_out(n_sweeps);
    std::vector<int64_t> cluster_out(n_sweeps);

    {
        py::gil_scoped_release release;

        seed_rng(seed_val);

        spin_tp s = spins.data();
        size_tp nlen_ptr = graph.nlen.data();
        NodesEdges ne = graph.node_edges.data();

        for (size_t step = 0; step < n_sweeps; ++step) {
            ene_out[step]  = calc_totEnergy(N, s, nlen_ptr, ne);
            magn_out[step] = calc_magn(N, s);
            cluster_out[step] = static_cast<int64_t>(
                sw_step(N, T, s, nlen_ptr, ne));
        }
    }

    py::array_t<int8_t> out_spins(N);
    std::memcpy(out_spins.mutable_data(), spins.data(), N * sizeof(int8_t));

    py::array_t<double> out_ene(n_sweeps);
    std::memcpy(out_ene.mutable_data(), ene_out.data(), n_sweeps * sizeof(double));

    py::array_t<double> out_magn(n_sweeps);
    std::memcpy(out_magn.mutable_data(), magn_out.data(), n_sweeps * sizeof(double));

    py::array_t<int64_t> out_clusters(n_sweeps);
    std::memcpy(out_clusters.mutable_data(), cluster_out.data(), n_sweeps * sizeof(int64_t));

    return py::make_tuple(out_spins, out_ene, out_magn, out_clusters);
}


/* ==================================================================
 * Topological Metropolis: MC in spectral coefficient space
 * ==================================================================
 *
 * Performs Metropolis updates by perturbing coefficients of a spectral
 * decomposition.  field = V^T @ c  is maintained incrementally; spins
 * are binarized from the field and energy is evaluated on the graph.
 */

static inline double rng_normal() {
    /* Box-Muller: two uniform → one normal deviate. */
    double u1 = RNG_dbl() + 1e-300;      /* avoid log(0) */
    double u2 = RNG_dbl();
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
}

static py::tuple topo_met_sampling(
    py::array_t<int8_t, py::array::c_style> spins_in,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr,
    const py::array_t<double, py::array::c_style>& subspace_vecs_in,   /* (M, N) */
    const py::array_t<double>& init_coefficients_in,                   /* (M,)   */
    const py::array_t<double>& init_field_in,                          /* (N,)   */
    const py::array_t<double>& mode_weights_in,                        /* (M,)   */
    double T,
    size_t n_steps,
    double sigma_init,
    size_t chunk_size,
    bool   do_polish,
    size_t polish_sweeps,
    uint64_t seed_val
) {
    /* ---- Parse inputs ------------------------------------------------ */
    auto s_buf = spins_in.request();
    size_t N = static_cast<size_t>(s_buf.size);

    auto sv_buf = subspace_vecs_in.request();
    if (sv_buf.ndim != 2)
        throw std::runtime_error("subspace_vecs must be 2-D (M, N)");
    size_t M = static_cast<size_t>(sv_buf.shape[0]);
    if (static_cast<size_t>(sv_buf.shape[1]) != N)
        throw std::runtime_error("subspace_vecs cols must equal N");

    auto coeff_buf  = init_coefficients_in.request();
    auto field_buf  = init_field_in.request();
    auto mw_buf     = mode_weights_in.request();

    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    /* ---- Copy mutable state ------------------------------------------ */
    std::vector<int8_t> spins(N);
    std::memcpy(spins.data(), s_buf.ptr, N * sizeof(int8_t));

    const double *V = static_cast<const double*>(sv_buf.ptr);   /* row-major (M,N) */
    std::vector<double> coeffs(M);
    std::memcpy(coeffs.data(), coeff_buf.ptr, M * sizeof(double));
    std::vector<double> field(N);
    std::memcpy(field.data(), field_buf.ptr, N * sizeof(double));
    std::vector<double> mw(M);
    std::memcpy(mw.data(), mw_buf.ptr, M * sizeof(double));

    /* ---- Pre-compute mode CDF for weighted selection ---- */
    std::vector<double> cdf(M);
    cdf[0] = mw[0];
    for (size_t k = 1; k < M; ++k) cdf[k] = cdf[k-1] + mw[k];

    /* ---- Adaptive sigma per mode ---- */
    std::vector<double> sigma(M, sigma_init);
    std::vector<size_t> accept_cnt(M, 0);
    std::vector<size_t> total_cnt(M, 0);

    /* ---- Output arrays ---- */
    std::vector<double> ene_out(n_steps);
    std::vector<double> magn_out(n_steps);

    /* ---- Best tracking ---- */
    std::vector<int8_t> best_spins(spins);
    double best_E = calc_totEnergy(N, spins.data(), graph.nlen.data(),
                                   graph.node_edges.data());
    double E_current = best_E;

    /* ---- Temporary buffers ---- */
    std::vector<double> field_new(N);
    std::vector<int8_t> spins_new(N);

    {
        py::gil_scoped_release release;
        seed_rng(seed_val);

        for (size_t step = 0; step < n_steps; ++step) {
            /* 1 MC step = N coefficient proposals */
            for (size_t prop = 0; prop < N; ++prop) {
                /* 1. Select mode via CDF */
                double r = RNG_dbl();
                size_t m = 0;
                for (size_t k = 0; k < M; ++k) {
                    if (r < cdf[k]) { m = k; break; }
                }

                /* 2. Propose delta ~ N(0, sigma[m]) */
                double delta = rng_normal() * sigma[m];

                /* 3. Incremental field update + binarize */
                const double *vm = V + m * N;
                for (size_t i = 0; i < N; ++i) {
                    field_new[i] = field[i] + delta * vm[i];
                    double f = field_new[i];
                    if      (f > 0.0)  spins_new[i] =  1;
                    else if (f < 0.0)  spins_new[i] = -1;
                    else               spins_new[i] = spins[i];
                }

                /* 4. Compute energy of proposed config */
                double E_new = calc_totEnergy(N, spins_new.data(),
                                              graph.nlen.data(),
                                              graph.node_edges.data());
                double dE = E_new - E_current;

                /* 5. Metropolis accept/reject */
                bool accept = (dE <= 0.0);
                if (!accept && T > 0.0) {
                    accept = (RNG_dbl() < std::exp(-dE / T));
                }

                if (accept) {
                    coeffs[m] += delta;
                    std::swap(field, field_new);
                    std::swap(spins, spins_new);
                    E_current = E_new;
                    accept_cnt[m]++;
                }

                total_cnt[m]++;
            }

            /* Record once per MC step (after N proposals) */
            ene_out[step]  = E_current;
            magn_out[step] = calc_magn(N, spins.data());

            if (E_current < best_E) {
                best_E = E_current;
                std::memcpy(best_spins.data(), spins.data(), N);
            }

            /* 6. Adaptive sigma every chunk_size steps */
            if (chunk_size > 0 && (step + 1) % chunk_size == 0) {
                for (size_t mm = 0; mm < M; ++mm) {
                    if (total_cnt[mm] > 0) {
                        double rate = (double)accept_cnt[mm] / total_cnt[mm];
                        if (rate > 0.35) sigma[mm] *= 1.1;
                        else if (rate < 0.25) sigma[mm] *= 0.9;
                    }
                }
                std::fill(accept_cnt.begin(), accept_cnt.end(), 0);
                std::fill(total_cnt.begin(), total_cnt.end(), 0);

                /* Optional T=0 greedy polish */
                if (do_polish) {
                    /* Use a zero-field array for polish */
                    std::vector<double> zero_h(N, 0.0);
                    initialize_glauberMetropolis(0.0, zero_h.data());
                    std::vector<int8_t> polish_s(spins);
                    for (size_t ps = 0; ps < polish_sweeps; ++ps) {
                        bool any_flip = false;
                        for (size_t i = 0; i < N; ++i) {
                            int8_t old_s = polish_s[i];
                            glauberMetropolis_1step_T0(
                                i, 0.0, polish_s.data(),
                                graph.nlen[i],
                                graph.node_edges[i]
                            );
                            if (polish_s[i] != old_s) any_flip = true;
                        }
                        if (!any_flip) break;
                    }
                    double E_pol = calc_totEnergy(
                        N, polish_s.data(),
                        graph.nlen.data(), graph.node_edges.data()
                    );
                    if (E_pol < best_E) {
                        best_E = E_pol;
                        best_spins = polish_s;
                    }
                }
            }
        }
    }

    /* ---- Convert to numpy ---- */
    auto out_s = py::array_t<int8_t>(N);
    std::memcpy(out_s.mutable_data(), spins.data(), N);
    auto out_best_s = py::array_t<int8_t>(N);
    std::memcpy(out_best_s.mutable_data(), best_spins.data(), N);
    auto out_coeffs = py::array_t<double>(M);
    std::memcpy(out_coeffs.mutable_data(), coeffs.data(), M * sizeof(double));

    auto out_ene  = py::array_t<double>(n_steps);
    auto out_magn = py::array_t<double>(n_steps);
    std::memcpy(out_ene.mutable_data(), ene_out.data(), n_steps * sizeof(double));
    std::memcpy(out_magn.mutable_data(), magn_out.data(), n_steps * sizeof(double));

    return py::make_tuple(out_s, out_ene, out_magn, out_best_s, best_E, out_coeffs);
}


/* ==================================================================
 * Cross-Entropy Method: population-based spectral optimizer
 * ==================================================================
 *
 * Samples K coefficient vectors c ~ N(mu, diag(sigma^2)),
 * maps to spins via sign(V @ c), optionally applies greedy T=0
 * quench, evaluates energy, selects elites, updates distribution.
 * Multi-restart keeps the global best.
 */

/* Per-site energy: E/N = -sum_{i<j} w_ij * s_i * s_j / N
 * Same convention as calc_totEnergy used by SA/Metropolis/PT. */
static inline double persite_edge_energy(
    size_t N,
    const int8_t *s,
    const size_t *nlen,
    const NodeEdges *ne
) {
    double E = 0.0;
    for (size_t i = 0; i < N; ++i) {
        double si = s[i];
        for (size_t p = 0; p < nlen[i]; ++p) {
            size_t j = ne[i].neighbors[p];
            if (j > i) {
                E -= ne[i].weights[p] * si * s[j];
            }
        }
    }
    return E / N;
}

/* T=0 greedy descent using graph CSR (no external field). */
static inline void greedy_quench(
    size_t N,
    int8_t *s,
    const size_t *nlen,
    const NodeEdges *ne,
    size_t max_sweeps
) {
    for (size_t sw = 0; sw < max_sweeps; ++sw) {
        bool any_flip = false;
        for (size_t i = 0; i < N; ++i) {
            double h_eff = 0.0;
            for (size_t p = 0; p < nlen[i]; ++p) {
                h_eff += ne[i].weights[p] * s[ne[i].neighbors[p]];
            }
            double dE = 2.0 * s[i] * h_eff;
            if (dE < 0.0) {
                s[i] = -s[i];
                any_flip = true;
            }
        }
        if (!any_flip) break;
    }
}

static py::tuple topo_cem_sampling(
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr,
    const py::array_t<double, py::array::c_style>& subspace_vecs_in,  /* (M, N) */
    size_t N,
    size_t cem_iter,
    size_t pop_size,
    double elite_frac,
    double init_sigma,
    double smoothing,
    double sigma_floor,
    double sigma_ceiling,
    size_t n_restarts,
    bool   do_greedy,
    size_t greedy_sweeps,
    bool   do_final_polish,
    size_t polish_sweeps,
    uint64_t seed_val
) {
    /* ---- Parse subspace ---- */
    auto sv_buf = subspace_vecs_in.request();
    if (sv_buf.ndim != 2)
        throw std::runtime_error("subspace_vecs must be 2-D (M, N)");
    size_t M = static_cast<size_t>(sv_buf.shape[0]);
    if (static_cast<size_t>(sv_buf.shape[1]) != N)
        throw std::runtime_error("subspace_vecs cols must equal N");
    const double *V = static_cast<const double*>(sv_buf.ptr);

    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    size_t n_elite = static_cast<size_t>(elite_frac * pop_size);
    if (n_elite < 1) n_elite = 1;

    /* ---- Output arrays ---- */
    std::vector<double> ene_trace(cem_iter);
    std::vector<double> restart_energies(n_restarts);

    double global_best_E = 1e30;
    std::vector<int8_t> global_best_spins(N, 1);
    std::vector<double> global_best_coeffs(M, 0.0);

    /* ---- Population buffers ---- */
    std::vector<double> coeffs_pop(pop_size * M);
    std::vector<int8_t> spins_pop(pop_size * N);
    std::vector<double> energies(pop_size);
    std::vector<size_t> sorted_idx(pop_size);

    {
        py::gil_scoped_release release;
        seed_rng(seed_val);

        for (size_t r = 0; r < n_restarts; ++r) {
            /* Initialize distribution */
            std::vector<double> mu(M, 0.0);
            std::vector<double> sigma(M, init_sigma);

            double restart_best_E = 1e30;
            std::vector<int8_t> restart_best_spins(N, 1);
            std::vector<double> restart_best_coeffs(M, 0.0);

            for (size_t it = 0; it < cem_iter; ++it) {
                /* 1. Sample K coefficient vectors: c ~ N(mu, diag(sigma^2)) */
                for (size_t k = 0; k < pop_size; ++k) {
                    double *ck = coeffs_pop.data() + k * M;
                    for (size_t j = 0; j < M; ++j) {
                        ck[j] = mu[j] + sigma[j] * rng_normal();
                    }
                }

                /* 2. Map to spins: sign(C @ V) and optionally greedy quench */
                for (size_t k = 0; k < pop_size; ++k) {
                    const double *ck = coeffs_pop.data() + k * M;
                    int8_t *sk = spins_pop.data() + k * N;

                    /* field_i = sum_j ck[j] * V[j*N + i] */
                    for (size_t i = 0; i < N; ++i) {
                        double f = 0.0;
                        for (size_t j = 0; j < M; ++j) {
                            f += ck[j] * V[j * N + i];
                        }
                        sk[i] = (f > 0.0) ? 1 : ((f < 0.0) ? -1 : 1);
                    }

                    /* Optional greedy quench */
                    if (do_greedy) {
                        greedy_quench(N, sk, graph.nlen.data(),
                                      graph.node_edges.data(), greedy_sweeps);
                    }

                    /* Evaluate energy */
                    energies[k] = persite_edge_energy(N, sk, graph.nlen.data(),
                                                  graph.node_edges.data());
                }

                /* 3. Sort by energy (ascending) */
                for (size_t k = 0; k < pop_size; ++k) sorted_idx[k] = k;
                std::sort(sorted_idx.begin(), sorted_idx.end(),
                          [&](size_t a, size_t b) { return energies[a] < energies[b]; });

                /* 4. Elite statistics */
                std::vector<double> new_mu(M, 0.0);
                std::vector<double> new_var(M, 0.0);
                for (size_t e = 0; e < n_elite; ++e) {
                    const double *ck = coeffs_pop.data() + sorted_idx[e] * M;
                    for (size_t j = 0; j < M; ++j) new_mu[j] += ck[j];
                }
                for (size_t j = 0; j < M; ++j) new_mu[j] /= n_elite;

                for (size_t e = 0; e < n_elite; ++e) {
                    const double *ck = coeffs_pop.data() + sorted_idx[e] * M;
                    for (size_t j = 0; j < M; ++j) {
                        double d = ck[j] - new_mu[j];
                        new_var[j] += d * d;
                    }
                }
                for (size_t j = 0; j < M; ++j) {
                    double new_sig = std::sqrt(new_var[j] / n_elite) + 1e-12;
                    mu[j] = smoothing * mu[j] + (1.0 - smoothing) * new_mu[j];
                    sigma[j] = smoothing * sigma[j] + (1.0 - smoothing) * new_sig;
                    if (sigma[j] < sigma_floor)   sigma[j] = sigma_floor;
                    if (sigma[j] > sigma_ceiling)  sigma[j] = sigma_ceiling;
                }

                /* Track best this iteration */
                size_t best_k = sorted_idx[0];
                if (energies[best_k] < restart_best_E) {
                    restart_best_E = energies[best_k];
                    std::memcpy(restart_best_spins.data(),
                                spins_pop.data() + best_k * N, N);
                    std::memcpy(restart_best_coeffs.data(),
                                coeffs_pop.data() + best_k * M, M * sizeof(double));
                }

                ene_trace[it] = restart_best_E;
            }

            /* Optional final polish */
            if (do_final_polish) {
                greedy_quench(N, restart_best_spins.data(),
                              graph.nlen.data(), graph.node_edges.data(),
                              polish_sweeps);
                double E_pol = persite_edge_energy(N, restart_best_spins.data(),
                                               graph.nlen.data(),
                                               graph.node_edges.data());
                if (E_pol < restart_best_E) restart_best_E = E_pol;
            }

            restart_energies[r] = restart_best_E;

            if (restart_best_E < global_best_E) {
                global_best_E = restart_best_E;
                global_best_spins = restart_best_spins;
                global_best_coeffs = restart_best_coeffs;
            }
        }
    }

    /* ---- Convert to numpy ---- */
    auto out_spins = py::array_t<int8_t>(N);
    std::memcpy(out_spins.mutable_data(), global_best_spins.data(), N);

    auto out_coeffs = py::array_t<double>(M);
    std::memcpy(out_coeffs.mutable_data(), global_best_coeffs.data(),
                M * sizeof(double));

    auto out_restart_E = py::array_t<double>(n_restarts);
    std::memcpy(out_restart_E.mutable_data(), restart_energies.data(),
                n_restarts * sizeof(double));

    auto out_ene = py::array_t<double>(cem_iter);
    std::memcpy(out_ene.mutable_data(), ene_trace.data(),
                cem_iter * sizeof(double));

    return py::make_tuple(out_spins, global_best_E, out_coeffs,
                          out_restart_E, out_ene);
}


static double calc_energy(
    py::array_t<int8_t, py::array::c_style> spins_in,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr
) {
    auto s_buf = spins_in.request();
    size_t N = static_cast<size_t>(s_buf.size);

    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    return calc_totEnergy(
        N,
        static_cast<spin_tp>(s_buf.ptr),
        graph.nlen.data(),
        graph.node_edges.data()
    );
}

/* ==================================================================
 * Module definition
 * ================================================================== */
PYBIND11_MODULE(_ising_native, m) {
    m.doc() = R"pbdoc(
        Native Ising dynamics kernels via pybind11.

        Wraps the C Glauber-Metropolis, Simulated Annealing, and
        Parallel Tempering kernels for direct Python invocation
        without subprocess or file I/O.
    )pbdoc";

    m.def("metropolis_sampling", &metropolis_sampling,
        R"pbdoc(
Run Metropolis sampling at fixed temperature.

Parameters
----------
spins : ndarray[int8]
    Initial spin configuration (N,).
neigh_indices : ndarray[int64]
    Concatenated neighbor indices.
neigh_weights : ndarray[float64]
    Concatenated neighbor edge weights (signs).
neigh_ptr : ndarray[int64]
    CSR row pointers (N+1,).
field : ndarray[float64]
    External field (N,).
T : float
    Temperature.
n_sweeps : int
    Number of MC sweeps.
seed : int
    RNG seed.
update_mode : str
    "sequential" or "asynchronous".

Returns
-------
tuple[ndarray, ndarray, ndarray]
    (final_spins, energy_trace, magnetization_trace)
        )pbdoc",
        py::arg("spins"), py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("field"), py::arg("T"),
        py::arg("n_sweeps"), py::arg("seed"),
        py::arg("update_mode") = "asynchronous");

    m.def("sa_sampling", &sa_sampling,
        R"pbdoc(
Run simulated annealing.

Returns
-------
tuple[ndarray, ndarray, ndarray, ndarray]
    (final_spins, energy_trace, magnetization_trace, temperature_schedule)
        )pbdoc",
        py::arg("spins"), py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("field"),
        py::arg("T_init"), py::arg("T_final"),
        py::arg("schedule"), py::arg("cooling_rate"),
        py::arg("steps_per_T"), py::arg("n_temps"),
        py::arg("seed"), py::arg("update_mode") = "asynchronous");

    m.def("pt_sampling", &pt_sampling,
        R"pbdoc(
Run parallel tempering (replica exchange Monte Carlo).

Returns
-------
tuple[ndarray, ndarray, ndarray, ndarray]
    (final_spins, energy[n_rep,n_ex], magn[n_rep,n_ex], exchanges[n_ex,n_rep-1])
        )pbdoc",
        py::arg("spins"), py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("field"),
        py::arg("T_ladder"), py::arg("steps_per_exchange"),
        py::arg("n_exchanges"), py::arg("seed"),
        py::arg("update_mode") = "asynchronous");

    m.def("wolff_sampling", &wolff_sampling,
        R"pbdoc(
Run Wolff cluster algorithm at fixed temperature.

Returns
-------
tuple[ndarray, ndarray, ndarray, ndarray]
    (final_spins, energy_trace, magnetization_trace, cluster_sizes)
        )pbdoc",
        py::arg("spins"), py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("field"), py::arg("T"),
        py::arg("n_sweeps"), py::arg("seed"));

    m.def("sw_sampling", &sw_sampling,
        R"pbdoc(
Run Swendsen-Wang cluster algorithm at fixed temperature.

Returns
-------
tuple[ndarray, ndarray, ndarray, ndarray]
    (final_spins, energy_trace, magnetization_trace, n_clusters_per_step)
        )pbdoc",
        py::arg("spins"), py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("field"), py::arg("T"),
        py::arg("n_sweeps"), py::arg("seed"));

    m.def("topo_met_sampling", &topo_met_sampling,
        R"pbdoc(
Run Topological Metropolis: MC in spectral coefficient space.

Parameters
----------
spins : ndarray[int8]
    Initial spins (N,).
neigh_indices, neigh_weights, neigh_ptr : CSR graph arrays.
subspace_vecs : ndarray[float64]
    Eigenvectors (M, N) row-major.
init_coefficients : ndarray[float64]
    Initial spectral coefficients (M,).
init_field : ndarray[float64]
    Initial continuous field (N,).
mode_weights : ndarray[float64]
    Mode selection probabilities (M,).
T : float
    Metropolis temperature.
n_steps : int
    Number of MC steps (coordinate-wise).
sigma_init : float
    Initial proposal std per mode.
chunk_size : int
    Adaptive sigma adjustment interval.
do_polish : bool
    Greedy T=0 polish after chunks.
polish_sweeps : int
    Max sweeps for polish.
seed : int
    RNG seed.

Returns
-------
tuple
    (final_spins, energy_trace, magn_trace,
     best_spins, best_energy, final_coefficients)
        )pbdoc",
        py::arg("spins"), py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("subspace_vecs"),
        py::arg("init_coefficients"), py::arg("init_field"),
        py::arg("mode_weights"),
        py::arg("T"), py::arg("n_steps"),
        py::arg("sigma_init"), py::arg("chunk_size"),
        py::arg("do_polish"), py::arg("polish_sweeps"),
        py::arg("seed"));

    m.def("topo_cem_sampling", &topo_cem_sampling,
        R"pbdoc(
Run Cross-Entropy Method spectral optimizer.

Parameters
----------
neigh_indices, neigh_weights, neigh_ptr : CSR graph arrays.
subspace_vecs : ndarray[float64]
    Eigenvectors (M, N) row-major.
N : int
    Number of nodes.
cem_iter : int
    CEM iterations per restart.
pop_size : int
    Population size K.
elite_frac : float
    Fraction of elites.
init_sigma : float
    Initial coefficient std.
smoothing : float
    Exponential smoothing for mu/sigma updates.
sigma_floor, sigma_ceiling : float
    Sigma bounds.
n_restarts : int
    Number of independent CEM restarts.
do_greedy : bool
    Greedy T=0 quench per sample.
greedy_sweeps : int
    Max sweeps per greedy quench.
do_final_polish : bool
    Final polish of restart best.
polish_sweeps : int
    Max sweeps for final polish.
seed : int
    RNG seed.

Returns
-------
tuple
    (best_spins, best_energy, best_coeffs, restart_energies, ene_trace)
        )pbdoc",
        py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("subspace_vecs"),
        py::arg("N"),
        py::arg("cem_iter"), py::arg("pop_size"),
        py::arg("elite_frac"), py::arg("init_sigma"),
        py::arg("smoothing"), py::arg("sigma_floor"),
        py::arg("sigma_ceiling"), py::arg("n_restarts"),
        py::arg("do_greedy"), py::arg("greedy_sweeps"),
        py::arg("do_final_polish"), py::arg("polish_sweeps"),
        py::arg("seed"));

    m.def("calc_energy", &calc_energy,
        R"pbdoc(
Calculate total Ising energy for a spin configuration.

Returns
-------
float
    Total energy E = -sum_{(i,j)} w_ij * s_i * s_j
        )pbdoc",
        py::arg("spins"), py::arg("neigh_indices"),
        py::arg("neigh_weights"), py::arg("neigh_ptr"));
}
