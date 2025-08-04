#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include <random>

namespace py = pybind11;

std::vector<int> random_walk(const py::array_t<double>& data,
                             const py::array_t<int>& indices,
                             const py::array_t<int>& indptr,
                             int num_steps, int start_node) {
    // Read the CSR matrix
    auto data_unchecked = data.unchecked<1>();
    auto indices_unchecked = indices.unchecked<1>();
    auto indptr_unchecked = indptr.unchecked<1>();

    std::vector<int> signs(num_steps);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, data_unchecked.size() - 1);

    int current_node = start_node;
    int current_sign = 1; // Start with a positive sign

    for (int step = 0; step < num_steps; ++step) {
        int row_start = indptr_unchecked(current_node);
        int row_end = indptr_unchecked(current_node + 1);

        if (row_start == row_end) {
            // No neighbors to walk to
            signs[step] = current_sign;
            continue;
        }

        int next_idx = distrib(gen) % (row_end - row_start) + row_start;
        int next_node = indices_unchecked(next_idx);
        double link_value = data_unchecked(next_idx);

        if (link_value < 0) {
            current_sign = -current_sign; // Switch sign if the link is negative
        }

        current_node = next_node;
        signs[step] = current_sign;
    }

    return signs;
}

PYBIND11_MODULE(random_walk, m) {
    m.def("random_walk", &random_walk, "Perform a random walk on the adjacency matrix",
          py::arg("data"), py::arg("indices"), py::arg("indptr"), py::arg("num_steps"), py::arg("start_node"));
}
