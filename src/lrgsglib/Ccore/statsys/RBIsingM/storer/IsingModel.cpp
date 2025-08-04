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
#include <random>

namespace py = boost::python;

#define ENABLE_TIMING

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

class IsingModel {
public:
    IsingModel(py::list edge_list, py::list sign_list, int width, int height, double temperature, std::string mode = "sequential")
        : width(width), height(height), temperature(temperature), mode(mode) {
        rng.seed(std::random_device{}());
#ifdef ENABLE_TIMING
        auto start = std::chrono::high_resolution_clock::now();
#endif
        initializeNeighborsAndSigns(edge_list, sign_list);
#ifdef ENABLE_TIMING
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "initializeNeighborsAndSigns took: " << elapsed.count() << " seconds\n";
#endif
        initializeSpins();
    }

    void simulate(int steps, int frame_rate) {
        int frame_interval = steps / frame_rate;
#ifdef ENABLE_TIMING
        auto start = std::chrono::high_resolution_clock::now();
#endif
        if (mode == "sequential") {
            for (int i = 0; i < steps; ++i) {
                sequentialUpdate();
                if (i % frame_interval == 0) {
                    captureFrame(i);
                }
            }
        } else if (mode == "synchronous") {
            for (int i = 0; i < steps; ++i) {
                synchronousUpdate();
                if (i % frame_interval == 0) {
                    captureFrame(i);
                }
            }
        } else if (mode == "asynchronous") {
            for (int i = 0; i < steps; ++i) {
                asynchronousUpdate();
                if (i % frame_interval == 0) {
                    captureFrame(i);
                }
            }
        } else {
            throw std::invalid_argument("Invalid mode. Use 'sequential', 'synchronous' or 'asynchronous'.");
        }
#ifdef ENABLE_TIMING
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "simulate took: " << elapsed.count() << " seconds\n";
#endif
    }

    py::list getSpins() const {
#ifdef ENABLE_TIMING
        auto start = std::chrono::high_resolution_clock::now();
#endif
        py::list spin_list;
        for (const auto& spin : spins) {
            spin_list.append(spin);
        }
#ifdef ENABLE_TIMING
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "getSpins took: " << elapsed.count() << " seconds\n";
#endif
        return spin_list;
    }

    double getEnergy() const {
        double energy = 0.0;
        for (int node = 0; node < width * height; ++node) {
            for (const auto& neighbor : neighbors.at(node)) {
                auto edge = normalize_edge(node, neighbor);
                energy += -spins[node] * spins[neighbor] * signs.at(edge);
            }
        }
        return energy / 2.0; // Each pair counted twice
    }

    double getMagnetization() const {
        double magnetization = 0.0;
        for (const auto& spin : spins) {
            magnetization += spin;
        }
        return magnetization / (width * height);
    }

    py::list getFrameSpins() const {
        return frame_spins;
    }

    py::list getFrameEnergies() const {
        return frame_energies;
    }

    py::list getFrameMagnetizations() const {
        return frame_magnetizations;
    }

private:
    std::map<std::pair<int, int>, int> signs;
    std::map<int, std::vector<int>> neighbors;
    std::vector<int> spins;
    int width;
    int height;
    double temperature;
    std::string mode;

    std::mt19937 rng;
    std::uniform_real_distribution<> dist{0.0, 1.0};

    py::list frame_spins;
    py::list frame_energies;
    py::list frame_magnetizations;

    void initializeNeighborsAndSigns(py::list edge_list, py::list sign_list) {
        std::vector<std::pair<int, int>> edges;
        std::vector<int> signs_vec;

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

    void initializeSpins() {
#ifdef ENABLE_TIMING
        auto start = std::chrono::high_resolution_clock::now();
#endif
        spins.resize(width * height, 1);
        for (auto& spin : spins) {
            spin = (dist(rng) > 0.5) ? 1 : -1; // Randomly set spin to +1 or -1
        }
#ifdef ENABLE_TIMING
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "initializeSpins took: " << elapsed.count() << " seconds\n";
#endif
    }

    double calculateEnergyChange(int node) const {
        double deltaE = 0.0;
        for (const auto& neighbor : neighbors.at(node)) {
            auto edge = normalize_edge(node, neighbor);
            deltaE += 2 * spins[node] * spins[neighbor] * signs.at(edge);
        }
        return deltaE;
    }

    void sequentialUpdate() {
#ifdef ENABLE_TIMING
        auto start = std::chrono::high_resolution_clock::now();
#endif
        for (int node = 0; node < width * height; ++node) {
            double deltaE = calculateEnergyChange(node);
            if (deltaE <= 0 || dist(rng) < std::exp(-deltaE / temperature)) {
                spins[node] *= -1;
            }
        }
#ifdef ENABLE_TIMING
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "sequentialUpdate took: " << elapsed.count() << " seconds\n";
#endif
    }

    void synchronousUpdate() {
#ifdef ENABLE_TIMING
        auto start = std::chrono::high_resolution_clock::now();
#endif
        std::vector<int> new_spins = spins;
        for (int node = 0; node < width * height; ++node) {
            double deltaE = calculateEnergyChange(node);
            if (deltaE <= 0 || dist(rng) < std::exp(-deltaE / temperature)) {
                new_spins[node] *= -1;
            }
        }
        spins = new_spins;
#ifdef ENABLE_TIMING
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "synchronousUpdate took: " << elapsed.count() << " seconds\n";
#endif
    }

    void asynchronousUpdate() {
#ifdef ENABLE_TIMING
        auto start = std::chrono::high_resolution_clock::now();
#endif
        for (int step = 0; step < width * height; ++step) { // Ensure the same number of updates as sequential
            int node = rng() % (width * height);
            double deltaE = calculateEnergyChange(node);
            if (deltaE <= 0 || dist(rng) < std::exp(-deltaE / temperature)) {
                spins[node] *= -1;
            }
        }
#ifdef ENABLE_TIMING
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout << "asynchronousUpdate took: " << elapsed.count() << " seconds\n";
#endif
    }

    void captureFrame(int step) {
        py::list frame;
        for (const auto& spin : spins) {
            frame.append(spin);
        }
        frame_spins.append(frame);
        frame_energies.append(getEnergy());
        frame_magnetizations.append(getMagnetization());
    }
};

BOOST_PYTHON_MODULE(ising_model_store) {
    py::to_python_converter<std::pair<int, int>, pair_to_python_tuple>();
    pair_from_python_tuple();

    py::class_<IsingModel>("IsingModel", py::init<py::list, py::list, int, int, double, py::optional<std::string>>())
        .def("simulate", &IsingModel::simulate)
        .def("getSpins", &IsingModel::getSpins)
        .def("getEnergy", &IsingModel::getEnergy)
        .def("getMagnetization", &IsingModel::getMagnetization)
        .def("getFrameSpins", &IsingModel::getFrameSpins)
        .def("getFrameEnergies", &IsingModel::getFrameEnergies)
        .def("getFrameMagnetizations", &IsingModel::getFrameMagnetizations);
}
