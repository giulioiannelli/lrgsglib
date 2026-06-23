/**
 * @file rd_native.cpp
 * @brief Pybind11 wrapper for the reaction-diffusion RK4 ODE C kernel.
 *
 * In-process reaction-diffusion dynamics on signed CSR graphs -- no file I/O,
 * works with both NX and GT graphs. Reuses the shared rd_rhs + rk4_step kernels
 * verbatim (RK4 via the shared ContDynSys integrator, with per-step clamping to
 * [0, inf) matching ReactionDiffusionSimulator0.c). The dynamics are
 * deterministic (fixed diffusion + reaction + initial field), so py<->pb agree
 * numerically (not merely statistically), up to RK4-implementation tolerance.
 *
 * Naming convention (Python side):
 *   runlang = "pb" / "pb_rd"  -> rd_sampling()
 */
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

extern "C" {
#include "LRGSG_customs.h"
#include "LRGSG_contdynsys.h"
#include "LRGSG_rd.h"
}

namespace py = pybind11;

/* Build NodeEdges adjacency from flat numpy CSR arrays (mirrors the helper in
 * kuramoto_native.cpp / voter_native.cpp / potts_native.cpp). */
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

/* ==================================================================
 * Reaction-diffusion sampling: n_steps RK4 integration steps at fixed dt,
 * each followed by a clamp to [0, inf) (matches both run_py's
 * _apply_constraints and ReactionDiffusionSimulator0.c). Deterministic.
 *
 * Records the mean concentration AFTER each step when save_mean is set,
 * giving n_steps records.
 *
 * Returns (final_u[float64], mean_concentration[float64]).
 * ================================================================== */
static py::tuple rd_sampling(
    py::array_t<double, py::array::c_style> u_in,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr,
    double D, double dt,
    const std::string& reaction,
    double r, double a,
    size_t n_steps,
    bool save_mean
) {
    auto u_buf = u_in.request();
    size_t N = static_cast<size_t>(u_buf.size);

    std::vector<double> u(N);
    std::memcpy(u.data(), u_buf.ptr, N * sizeof(double));

    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    std::vector<double> mean_out;
    if (save_mean) mean_out.reserve(n_steps);

    {
        py::gil_scoped_release release;

        RDParams par;
        par.N          = N;
        par.D          = D;
        par.reaction   = parse_reaction_type(reaction.c_str());
        par.r          = r;
        par.a          = a;
        par.neigh_len  = graph.nlen.data();
        par.node_edges = graph.node_edges.data();

        double t = 0.0;
        for (size_t step = 0; step < n_steps; ++step) {
            rk4_step(N, u.data(), rd_rhs, t, dt, &par);
            for (size_t i = 0; i < N; ++i)
                if (u[i] < 0.0) u[i] = 0.0;
            t += dt;
            if (save_mean)
                mean_out.push_back(calc_mean_concentration(N, u.data()));
        }
    }

    py::array_t<double> out_u(N);
    std::memcpy(out_u.mutable_data(), u.data(), N * sizeof(double));

    py::array_t<double> out_mean(mean_out.size());
    if (!mean_out.empty())
        std::memcpy(out_mean.mutable_data(), mean_out.data(),
                    mean_out.size() * sizeof(double));

    return py::make_tuple(out_u, out_mean);
}

PYBIND11_MODULE(_rd_native, m) {
    m.doc() = "Native (pybind11) reaction-diffusion RK4 kernel on signed CSR graphs.";

    m.def("rd_sampling", &rd_sampling,
        R"pbdoc(
Integrate reaction-diffusion dynamics on a signed CSR graph (RK4).

Parameters
----------
u, neigh_indices, neigh_weights, neigh_ptr : ndarray
    Initial float64 concentration field and the CSR adjacency (symmetric,
    signed weights).
D, dt : diffusion constant and integration time step.
reaction : reaction type ("fisher_kpp", "bistable", "linear", "none").
r, a : reaction parameters (growth rate / bistable threshold).
n_steps : number of RK4 steps.
save_mean : record the mean concentration after each (clamped) step.

Returns
-------
tuple[ndarray, ndarray]
    (final_field[float64], mean_concentration[float64]); the mean array has
    length n_steps when save_mean, else 0.
        )pbdoc",
        py::arg("u"), py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("D"), py::arg("dt"),
        py::arg("reaction"), py::arg("r"), py::arg("a"),
        py::arg("n_steps"), py::arg("save_mean") = false);
}
