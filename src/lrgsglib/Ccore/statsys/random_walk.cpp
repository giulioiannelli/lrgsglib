#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <optional>
#include <random>
#include <vector>
#include <pybind11/stl.h>

namespace py = pybind11;

std::vector<int> signed_random_walk(
    const py::array_t<double, py::array::c_style | py::array::forcecast>& data,
    const py::array_t<int, py::array::c_style | py::array::forcecast>& indices,
    const py::array_t<int, py::array::c_style | py::array::forcecast>& indptr,
    int num_steps, int start_node, std::optional<std::uint64_t> seed) {
    if (num_steps < 0) {
        throw py::value_error("num_steps must be non-negative");
    }
    if (data.ndim() != 1 || indices.ndim() != 1 || indptr.ndim() != 1) {
        throw py::value_error("data, indices, and indptr must be 1D arrays");
    }
    if (data.size() != indices.size()) {
        throw py::value_error("data and indices must have the same length");
    }
    if (indptr.size() < 2) {
        throw py::value_error("indptr must have at least two entries");
    }

    const auto n_nodes = indptr.size() - 1;
    if (start_node < 0 || static_cast<py::ssize_t>(start_node) >= n_nodes) {
        throw py::value_error("start_node is out of bounds");
    }

    auto data_unchecked = data.unchecked<1>();
    auto indices_unchecked = indices.unchecked<1>();
    auto indptr_unchecked = indptr.unchecked<1>();

    std::mt19937 gen(seed.value_or(std::random_device{}()));
    std::vector<int> signs(static_cast<std::size_t>(num_steps));

    int current_node = start_node;
    int current_sign = 1;  // Start with a positive sign

    py::gil_scoped_release release;  // Heavy loop runs without the GIL
    for (int step = 0; step < num_steps; ++step) {
        const auto row_start = indptr_unchecked(current_node);
        const auto row_end = indptr_unchecked(current_node + 1);

        if (row_start == row_end) {
            // No neighbors to walk to
            signs[static_cast<std::size_t>(step)] = current_sign;
            continue;
        }

        std::uniform_int_distribution<py::ssize_t> neighbor_dist(
            row_start, row_end - 1);
        const auto next_idx = neighbor_dist(gen);
        const auto next_node = indices_unchecked(next_idx);
        const auto link_value = data_unchecked(next_idx);

        if (link_value < 0) {
            current_sign = -current_sign;  // Switch sign if the link is negative
        }

        current_node = next_node;
        signs[static_cast<std::size_t>(step)] = current_sign;
    }

    return signs;
}

PYBIND11_MODULE(_random_walk, m) {
    m.doc() = "C++ random walk on a signed CSR adjacency matrix";
    m.def(
        "random_walk", &signed_random_walk,
        R"pbdoc(
Perform a signed random walk on a CSR adjacency matrix.

Parameters
----------
data, indices, indptr: 1D numpy arrays
    CSR components of the adjacency matrix (weights can be positive or negative).
num_steps: int
    Number of steps to take.
start_node: int, default=0
    Starting node index.
seed: Optional[int], default=None
    RNG seed for reproducibility. If None, uses std::random_device.
        )pbdoc",
        py::arg("data"), py::arg("indices"), py::arg("indptr"),
        py::arg("num_steps"), py::arg("start_node") = 0,
        py::arg("seed") = py::none());
}
