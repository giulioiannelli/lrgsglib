#include <boost/python.hpp>
#include "pair_conversion.hpp"
#include "ising_model_store.hpp"

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
