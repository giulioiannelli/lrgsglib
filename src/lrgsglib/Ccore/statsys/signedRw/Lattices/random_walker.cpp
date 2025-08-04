#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <vector>
#include <map>
#include <cstdlib>
#include <ctime>
#include <utility>
#include <cmath>
#include <chrono>
#include <iostream>

#define ENABLE_TIMING
// #undef ENABLE_TIMING

namespace py = boost::python;

// Converter from Python tuple to std::pair<int, int>
struct pair_to_python_tuple {
    static PyObject* convert(const std::pair<int, int>& p) {
        return py::incref(py::make_tuple(p.first, p.second).ptr());
    }
};

struct pair_from_python_tuple {
    pair_from_python_tuple() {
        py::converter::registry::push_back(
            &convertible,
            &construct,
            py::type_id<std::pair<int, int>>());
    }

    static void* convertible(PyObject* obj_ptr) {
        if (!PyTuple_Check(obj_ptr)) return 0;
        if (PyTuple_Size(obj_ptr) != 2) return 0;
        return obj_ptr;
    }

    static void construct(PyObject* obj_ptr, py::converter::rvalue_from_python_stage1_data* data) {
        void* storage = ((py::converter::rvalue_from_python_storage<std::pair<int, int>>*)data)->storage.bytes;
        PyObject* first_obj = PyTuple_GetItem(obj_ptr, 0);
        PyObject* second_obj = PyTuple_GetItem(obj_ptr, 1);
        
        if (!PyLong_Check(first_obj) || !PyLong_Check(second_obj)) {
            PyErr_SetString(PyExc_TypeError, "Tuple elements must be integers");
            py::throw_error_already_set();
        }
        
        int first = static_cast<int>(PyLong_AsLong(first_obj));
        int second = static_cast<int>(PyLong_AsLong(second_obj));

        new (storage) std::pair<int, int>(first, second);
        data->convertible = storage;
    }
};

std::pair<int, int> normalize_edge(int a, int b) {
    return std::make_pair(std::min(a, b), std::max(a, b));
}

class RandomWalker {
public:
    RandomWalker(py::list edge_list, py::list sign_list, int width, int height, unsigned int coordination_number)
        : width(width), height(height), coordination_number(coordination_number) {
#ifdef ENABLE_TIMING
        auto start = std::chrono::high_resolution_clock::now();
#endif
        initializeNeighborsAndSigns(edge_list, sign_list);
#ifdef ENABLE_TIMING
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "initializeNeighborsAndSigns took: " << elapsed.count() << " seconds\n";
#endif
    }

    void walk(int start) {
#ifdef ENABLE_TIMING
        auto start_time = std::chrono::high_resolution_clock::now();
#endif

        int current = start;
        int start_x = current / width;
        int start_y = current % width;
        int current_x, current_y;
        distance = 0;
        euclidean_distance = 0.0;
        stopped_node = current; // Initialize the stopping node

        while (true) {
            const auto& neighbors = this->neighbors[current];
            if (isBorder(current)) break;

            int next = neighbors[rand() % neighbors.size()];
            auto normalized_edge = normalize_edge(current, next);
            if (signs[normalized_edge] == -1) break;

            current = next;
            stopped_node = current; // Update the stopping node
            distance++;
        }
        current_x = current / width;
        current_y = current % width;
        euclidean_distance = std::sqrt(std::pow(current_x - start_x, 2) + std::pow(current_y - start_y, 2));

#ifdef ENABLE_TIMING
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end_time - start_time;
        std::cout << "walk took: " << elapsed.count() << " seconds\n";
#endif
    }

    int get_step_distance() const {
        return distance;
    }

    double get_euclidean_distance() const {
        return euclidean_distance;
    }

    int get_stopped_node() const {
        return stopped_node;
    }

private:
    std::map<std::pair<int, int>, int> signs;
    std::map<int, std::vector<int>> neighbors;
    int distance = 0;
    double euclidean_distance = 0.0;
    int width;
    int height;
    unsigned int coordination_number;
    int stopped_node;

    void initializeNeighborsAndSigns(py::list edge_list, py::list sign_list) {
        std::vector<std::pair<int, int>> edges;
        std::vector<int> signs_vec;

        // Using Boost.Python's STL iterator to convert Python list to C++ vector
        std::copy(py::stl_input_iterator<std::pair<int, int>>(edge_list), py::stl_input_iterator<std::pair<int, int>>(), std::back_inserter(edges));
        std::copy(py::stl_input_iterator<int>(sign_list), py::stl_input_iterator<int>(), std::back_inserter(signs_vec));

        for (size_t i = 0; i < edges.size(); ++i) {
            auto edge = edges[i];
            int sign = signs_vec[i];
            auto normalized_edge = normalize_edge(edge.first, edge.second);
            neighbors[normalized_edge.first].push_back(normalized_edge.second);
            neighbors[normalized_edge.second].push_back(normalized_edge.first);
            this->signs[normalized_edge] = sign;
        }
    }

    bool isBorder(int node) {
        return neighbors[node].size() < coordination_number;
    }
};

BOOST_PYTHON_MODULE(random_walker) {
    py::to_python_converter<std::pair<int, int>, pair_to_python_tuple>();
    pair_from_python_tuple();

    py::class_<RandomWalker>("RandomWalker", py::init<py::list, py::list, int, int, unsigned int>())
        .def("walk", &RandomWalker::walk)
        .def("get_step_distance", &RandomWalker::get_step_distance)
        .def("get_euclidean_distance", &RandomWalker::get_euclidean_distance)
        .def("get_stopped_node", &RandomWalker::get_stopped_node);
}
