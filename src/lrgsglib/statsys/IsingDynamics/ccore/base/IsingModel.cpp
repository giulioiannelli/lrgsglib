#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <map>
#include <cstdlib>
#include <ctime>
#include <utility>
#include <cmath>
#include <chrono>
#include <iostream>
#include <random>

namespace py = pybind11;

// #define ENABLE_TIMING

static std::pair<int, int> normalize_edge(int a, int b) {
    return std::make_pair(std::min(a, b), std::max(a, b));
}

class IsingModel {
public:
    IsingModel(py::list edge_list, py::list sign_list, int width, int height, double temperature, std::string mode = "synchronous")
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

    void simulate(int steps) {
#ifdef ENABLE_TIMING
        auto start = std::chrono::high_resolution_clock::now();
#endif
        if (mode == "synchronous") {
            for (int i = 0; i < steps; ++i) {
                synchronousUpdate();
            }
        } else if (mode == "asynchronous") {
            for (int i = 0; i < steps; ++i) {
                asynchronousUpdate();
            }
        } else {
            throw std::invalid_argument("Invalid mode. Use 'synchronous' or 'asynchronous'.");
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

    void initializeNeighborsAndSigns(py::list edge_list, py::list sign_list) {
        for (size_t i = 0; i < static_cast<size_t>(py::len(edge_list)); ++i) {
            auto edge = edge_list[i].cast<std::pair<int, int>>();
            int sign = sign_list[i].cast<int>();
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

    void synchronousUpdate() {
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
        std::cout << "synchronousUpdate took: " << elapsed.count() << " seconds\n";
#endif
    }

    void asynchronousUpdate() {
#ifdef ENABLE_TIMING
        auto start = std::chrono::high_resolution_clock::now();
#endif
        for (int step = 0; step < width * height; ++step) { // Ensure the same number of updates as synchronous
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
};

PYBIND11_MODULE(ising_model, m) {
    m.doc() = "Ising model simulation with signed edges";

    py::class_<IsingModel>(m, "IsingModel")
        .def(py::init<py::list, py::list, int, int, double, std::string>(),
             py::arg("edge_list"), py::arg("sign_list"),
             py::arg("width"), py::arg("height"),
             py::arg("temperature"), py::arg("mode") = "synchronous")
        .def("simulate", &IsingModel::simulate)
        .def("getSpins", &IsingModel::getSpins);
}
