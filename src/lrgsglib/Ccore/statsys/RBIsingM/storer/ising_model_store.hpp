#ifndef ISING_MODEL_HPP
#define ISING_MODEL_HPP

#include <boost/python.hpp>
#include <vector>
#include <map>
#include <random>
#include <string>


namespace py = boost::python;

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

    py::list frame_spins;
    py::list frame_energies;
    py::list frame_magnetizations;
};

#endif // ISING_MODEL_HPP
