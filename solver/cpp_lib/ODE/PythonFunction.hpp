#ifndef PYTHONFUNCTION_HPP
#define PYTHONFUNCTION_HPP

#include <Python.h>
#include <numpy/ndarrayobject.h>
#include "Function.hpp"
#include "../utils.hpp"

class PythonFunction: public Function
{
public:
    PythonFunction(PyObject* callable): callable_(callable) {}

    virtual ~PythonFunction() {}

    virtual Real operator()(std::vector<Real>& arr, Real current_time) const
    {
        npy_intp dims[] = {arr.size()};
        PyObject* _arr = PyArray_New(&PyArray_Type, 1,
            dims, NPY_DOUBLE, NULL, &arr[0], 0, NPY_CARRAY, NULL);
        PyObject* _current_time = PyFloat_FromDouble(current_time);
        PyObject* _result = PyObject_CallFunctionObjArgs(callable_,
            _arr, _current_time, NULL);
        if (!_result)
        {
            // in case of error
            Py_DECREF(_current_time);
            Py_DECREF(_arr);
            translate_python_exception();
        }

        const Real retval = PyFloat_AsDouble(_result);

        Py_DECREF(_result);
        Py_DECREF(_current_time);
        Py_DECREF(_arr);
        return retval;
    }

private:
    PyObject* callable_;
};

#endif /* PYTHONFUNCTION_HPP */
