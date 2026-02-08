#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "ising_model_store.hpp"

namespace py = pybind11;

PYBIND11_MODULE(ising_model_store, m) {
    m.doc() = "Ising model simulation with frame storage for animation";

    py::class_<IsingModel>(m, "IsingModel")
        .def(py::init<py::list, py::list, int, int, double, std::string>(),
             py::arg("edge_list"), py::arg("sign_list"),
             py::arg("width"), py::arg("height"),
             py::arg("temperature"), py::arg("mode") = "sequential")
        .def("simulate", &IsingModel::simulate)
        .def("getSpins", &IsingModel::getSpins)
        .def("getEnergy", &IsingModel::getEnergy)
        .def("getMagnetization", &IsingModel::getMagnetization)
        .def("getFrameSpins", &IsingModel::getFrameSpins)
        .def("getFrameEnergies", &IsingModel::getFrameEnergies)
        .def("getFrameMagnetizations", &IsingModel::getFrameMagnetizations);
}
