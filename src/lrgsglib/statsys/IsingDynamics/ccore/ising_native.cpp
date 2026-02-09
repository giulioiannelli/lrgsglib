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
