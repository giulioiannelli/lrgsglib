#ifndef PAIR_CONVERSION_HPP
#define PAIR_CONVERSION_HPP

#include <boost/python.hpp>

namespace py = boost::python;

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

#endif // PAIR_CONVERSION_HPP
