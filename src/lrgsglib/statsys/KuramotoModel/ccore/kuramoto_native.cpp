/**
 * @file kuramoto_native.cpp
 * @brief Pybind11 wrapper for the Kuramoto RK4 ODE C kernel.
 *
 * In-process Kuramoto phase-oscillator dynamics on signed CSR graphs -- no file
 * I/O, works with both NX and GT graphs. Reuses the shared kuramoto_step_rk4 /
 * kuramoto_rhs kernels verbatim (RK4 + wrap-to-2pi via the shared ContDynSys
 * integrator). The dynamics are deterministic (fixed omega + initial phases),
 * so py<->pb agree numerically up to the RK4-implementation tolerance, not just
 * statistically.
 *
 * Naming convention (Python side):
 *   runlang = "pb" / "pb_kuramoto"  -> kuramoto_sampling()
 */
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

extern "C" {
#include "LRGSG_customs.h"
#include "LRGSG_contdynsys.h"
#include "LRGSG_kuramoto.h"
}

namespace py = pybind11;

/* Build NodeEdges adjacency from flat numpy CSR arrays (mirrors the helper in
 * voter_native.cpp / potts_native.cpp). */
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

/* Kuramoto order parameter r = |1/N sum exp(i theta)|; matches
 * KuramotoModel.order_parameter. */
static double order_parameter(size_t N, const double* theta) {
    double cs = 0.0, sn = 0.0;
    for (size_t i = 0; i < N; ++i) { cs += std::cos(theta[i]); sn += std::sin(theta[i]); }
    return std::sqrt(cs * cs + sn * sn) / static_cast<double>(N);
}

/* ==================================================================
 * Kuramoto sampling: n_steps RK4 integration steps at fixed dt.
 *
 * Records the order parameter AFTER each step (matching KuramotoModel.run_py)
 * when save_order_param is set, giving n_steps records. Deterministic.
 *
 * Returns (final_theta[float64], order_params[float64]).
 * ================================================================== */
static py::tuple kuramoto_sampling(
    py::array_t<double, py::array::c_style> theta_in,
    const py::array_t<int64_t>& neigh_indices,
    const py::array_t<double>&  neigh_weights,
    const py::array_t<int64_t>& neigh_ptr,
    py::array_t<double, py::array::c_style> omega_in,
    double K, double dt,
    size_t n_steps,
    bool save_order_param
) {
    auto th_buf = theta_in.request();
    size_t N = static_cast<size_t>(th_buf.size);

    std::vector<double> theta(N);
    std::memcpy(theta.data(), th_buf.ptr, N * sizeof(double));

    std::vector<double> omega(N);
    std::memcpy(omega.data(), omega_in.request().ptr, N * sizeof(double));

    GraphCSR graph(N, neigh_indices, neigh_weights, neigh_ptr);

    std::vector<double> order_out;
    if (save_order_param) order_out.reserve(n_steps);

    {
        py::gil_scoped_release release;

        KuramotoParams par;
        par.N = N;
        par.K = K;
        par.omega = omega.data();
        par.neigh_len = graph.nlen.data();
        par.node_edges = graph.node_edges.data();

        double t = 0.0;
        for (size_t step = 0; step < n_steps; ++step) {
            kuramoto_step_rk4(N, theta.data(), t, dt, &par);
            t += dt;
            if (save_order_param)
                order_out.push_back(order_parameter(N, theta.data()));
        }
    }

    py::array_t<double> out_theta(N);
    std::memcpy(out_theta.mutable_data(), theta.data(), N * sizeof(double));

    py::array_t<double> out_order(order_out.size());
    if (!order_out.empty())
        std::memcpy(out_order.mutable_data(), order_out.data(),
                    order_out.size() * sizeof(double));

    return py::make_tuple(out_theta, out_order);
}

PYBIND11_MODULE(_kuramoto_native, m) {
    m.doc() = "Native (pybind11) Kuramoto RK4 kernel on signed CSR graphs.";

    m.def("kuramoto_sampling", &kuramoto_sampling,
        R"pbdoc(
Integrate Kuramoto phase-oscillator dynamics on a signed CSR graph (RK4).

Parameters
----------
theta, neigh_indices, neigh_weights, neigh_ptr : ndarray
    Initial float64 phases and the CSR adjacency (symmetric, signed weights).
omega : ndarray
    Per-node natural frequencies (float64, length N).
K, dt : coupling strength and integration time step.
n_steps : number of RK4 steps.
save_order_param : record the order parameter r=|<e^{i theta}>| after each step.

Returns
-------
tuple[ndarray, ndarray]
    (final_phases[float64], order_parameter[float64]); the order array has
    length n_steps when save_order_param, else 0.
        )pbdoc",
        py::arg("theta"), py::arg("neigh_indices"), py::arg("neigh_weights"),
        py::arg("neigh_ptr"), py::arg("omega"), py::arg("K"), py::arg("dt"),
        py::arg("n_steps"), py::arg("save_order_param") = false);
}
