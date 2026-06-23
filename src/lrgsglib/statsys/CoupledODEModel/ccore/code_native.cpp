/**
 * @file code_native.cpp
 * @brief Pybind11 wrapper for the CoupledODEModel RK4 ODE C kernel.
 *
 * In-process generic coupled-ODE dynamics on signed CSR graphs -- no file I/O,
 * works with both NX and GT graphs. Reuses the shared code_rhs / rk4_step
 * kernels verbatim (RK4 via the shared ContDynSys integrator). The dynamics are
 * deterministic (fixed initial condition + RHS), so py<->pb agree numerically
 * up to the RK4-implementation tolerance, not just statistically.
 *
 * Naming convention (Python side):
 *   runlang = "pb" / "pb_code"  -> code_sampling()
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
#include "LRGSG_code.h"
}

namespace py = pybind11;

/* Build NodeEdges adjacency from flat numpy CSR arrays (mirrors the helper in
 * kuramoto_native.cpp / voter_native.cpp). */
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
 * Coupled-ODE sampling: n_steps RK4 integration steps at fixed dt.
 *
 * Wraps the existing code_rhs kernel verbatim through the shared rk4_step
 * integrator. Deterministic.
 *
 * Returns final_state[float64] (length N).
 * ================================================================== */
static py::array_t<double> code_sampling(
    py::array_t<double, py::array::c_style> x_in,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr,
    double coupling_strength,
    const std::string& coupling_type,
    const std::string& local_type,
    double a, double r, double K,
    double dt,
    size_t n_steps
) {
    auto x_buf = x_in.request();
    size_t N = static_cast<size_t>(x_buf.size);

    std::vector<double> x(N);
    std::memcpy(x.data(), x_buf.ptr, N * sizeof(double));

    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    {
        py::gil_scoped_release release;

        CODEParams par;
        par.N                 = N;
        par.coupling_strength = coupling_strength;
        par.coupling          = parse_coupling_type(coupling_type.c_str());
        par.local             = parse_local_type(local_type.c_str());
        par.a                 = a;
        par.r                 = r;
        par.K                 = K;
        par.neigh_len         = graph.nlen.data();
        par.node_edges        = graph.node_edges.data();

        double t = 0.0;
        for (size_t step = 0; step < n_steps; ++step) {
            rk4_step(N, x.data(), code_rhs, t, dt, &par);
            t += dt;
        }
    }

    py::array_t<double> out_x(N);
    std::memcpy(out_x.mutable_data(), x.data(), N * sizeof(double));
    return out_x;
}

PYBIND11_MODULE(_code_native, m) {
    m.doc() = "Native (pybind11) CoupledODEModel RK4 kernel on signed CSR graphs.";

    m.def("code_sampling", &code_sampling,
        R"pbdoc(
Integrate generic coupled-ODE dynamics on a signed CSR graph (RK4).

Parameters
----------
x, neigh_indices, neigh_weights, neigh_ptr : ndarray
    Initial float64 state and the CSR adjacency (symmetric, signed weights).
coupling_strength : float
    Multiplier for the coupling term.
coupling_type : {"linear"/"diffusive", "product"}
local_type : {"none", "linear", "logistic", "fitzhugh"/"fhn"}
a, r, K : local-function parameters (a: linear coeff, r: growth, K: capacity).
dt : integration time step.
n_steps : number of RK4 steps.

Returns
-------
ndarray
    final_state[float64] (length N).
        )pbdoc",
        py::arg("x"), py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("coupling_strength"),
        py::arg("coupling_type"), py::arg("local_type"),
        py::arg("a"), py::arg("r"), py::arg("K"),
        py::arg("dt"), py::arg("n_steps"));
}
