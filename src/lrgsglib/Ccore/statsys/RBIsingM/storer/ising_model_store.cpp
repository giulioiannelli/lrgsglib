#include "ising_model_store.hpp"
#include <chrono>
#include <cmath>
#include <iostream>
#include <algorithm>

std::pair<int, int> normalize_edge(int a, int b) {
    return std::make_pair(std::min(a, b), std::max(a, b));
}

IsingModel::IsingModel(py::list edge_list, py::list sign_list, int width, int height, double temperature, std::string mode)
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

void IsingModel::simulate(int steps, int frame_rate) {
#ifdef ENABLE_TIMING
    auto start = std::chrono::high_resolution_clock::now();
#endif
    if (mode == "sequential") {
        for (int i = 0; i < steps; ++i) {
            sequentialUpdate();
            if (i % frame_rate == 0) {
                captureFrame(i);
            }
        }
    } else if (mode == "synchronous") {
        for (int i = 0; i < steps; ++i) {
            synchronousUpdate();
            if (i % frame_rate == 0) {
                captureFrame(i);
            }
        }
    } else if (mode == "asynchronous") {
        for (int i = 0; i < steps; ++i) {
            asynchronousUpdate();
            if (i % frame_rate == 0) {
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

py::list IsingModel::getSpins() const {
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

double IsingModel::getEnergy() const {
    double energy = 0.0;
    for (int node = 0; node < width * height; ++node) {
        for (const auto& neighbor : neighbors.at(node)) {
            auto edge = normalize_edge(node, neighbor);
            energy += -spins[node] * spins[neighbor] * signs.at(edge);
        }
    }
    return energy / 2.0; // Each pair counted twice
}

double IsingModel::getMagnetization() const {
    double magnetization = 0.0;
    for (const auto& spin : spins) {
        magnetization += spin;
    }
    return magnetization / (width * height);
}

py::list IsingModel::getFrameSpins() const {
    return frame_spins;
}

py::list IsingModel::getFrameEnergies() const {
    return frame_energies;
}

py::list IsingModel::getFrameMagnetizations() const {
    return frame_magnetizations;
}

void IsingModel::initializeNeighborsAndSigns(py::list edge_list, py::list sign_list) {
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

    node_dist = std::uniform_int_distribution<>(0, width * height - 1);
}

void IsingModel::initializeSpins() {
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

double IsingModel::calculateEnergyChange(int node) const {
    double deltaE = 0.0;
    for (const auto& neighbor : neighbors.at(node)) {
        auto edge = normalize_edge(node, neighbor);
        deltaE += 2 * spins[node] * spins[neighbor] * signs.at(edge);
    }
    return deltaE;
}

void IsingModel::sequentialUpdate() {
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

void IsingModel::synchronousUpdate() {
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

void IsingModel::asynchronousUpdate() {
#ifdef ENABLE_TIMING
    auto start = std::chrono::high_resolution_clock::now();
#endif
    for (int step = 0; step < width * height; ++step) { // Ensure the same number of updates as sequential
        int node = node_dist(rng); // Use pre-defined distribution
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

void IsingModel::captureFrame(int step) {
    py::list frame;
    for (const auto& spin : spins) {
        frame.append(spin);
    }
    frame_spins.append(frame);
    frame_energies.append(getEnergy());
    frame_magnetizations.append(getMagnetization());
}
