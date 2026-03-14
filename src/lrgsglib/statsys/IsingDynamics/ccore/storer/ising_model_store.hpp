#ifndef ISING_MODEL_HPP
#define ISING_MODEL_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <map>
#include <random>
#include <string>


namespace py = pybind11;

std::pair<int, int> normalize_edge(int a, int b);

class IsingModel {
public:
    IsingModel(py::list edge_list, py::list sign_list, int width, int height, double temperature, std::string mode = "sequential");

    void simulate(int steps, int frame_rate);
    py::list getSpins() const;
    double getEnergy() const;
    double getMagnetization() const;
    py::list getFrameSpins() const;
    py::list getFrameEnergies() const;
    py::list getFrameMagnetizations() const;

private:
    void initializeNeighborsAndSigns(py::list edge_list, py::list sign_list);
    void initializeSpins();
    double calculateEnergyChange(int node) const;
    void sequentialUpdate();
    void synchronousUpdate();
    void asynchronousUpdate();
    void captureFrame(int step);

    std::map<std::pair<int, int>, int> signs;
    std::map<int, std::vector<int>> neighbors;
    std::vector<int> spins;
    int width;
    int height;
    double temperature;
    std::string mode;

    std::mt19937 rng;
    std::uniform_real_distribution<> dist{0.0, 1.0};
    std::uniform_int_distribution<> node_dist;

    std::vector<std::vector<int>> frame_spins;
    std::vector<double> frame_energies;
    std::vector<double> frame_magnetizations;
};

#endif // ISING_MODEL_HPP
