/*
 * for DAESolver
 */
#ifndef PYTHONFUNCTION_HPP
#define PYTHONFUNCTION_HPP

#include <Python.h>
#include <numpy/ndarrayobject.h>
#include "Function.hpp"
#include "../utils.hpp"

class PythonFunction
  :
  public Function
{
public:
    PythonFunction(PyObject* callable): callable_(callable) {}

    virtual ~PythonFunction() {}

    virtual Real operator()(std::vector<Real>& array_differential,
        std::vector<Real>& array_algebraic, Real current_time) const
    {
        npy_intp dims_differential[] = {array_differential.size()};
        npy_intp dims_algebraic[] = {array_algebraic.size()};
        PyObject* _array_differential = PyArray_New(&PyArray_Type, 1,
            dims_differential, NPY_DOUBLE, NULL, &array_differential[0],
            0, NPY_CARRAY, NULL);
        PyObject* _array_algebraic = PyArray_New(&PyArray_Type, 1,
            dims_algebraic, NPY_DOUBLE, NULL, &array_algebraic[0],
            0, NPY_CARRAY, NULL);
        PyObject* _current_time = PyFloat_FromDouble(current_time);
        PyObject* _result = PyObject_CallFunctionObjArgs(callable_,
            _array_differential, _array_algebraic, _current_time, NULL);
        if (!_result)
        {
            // in case of error
            Py_DECREF(_current_time);
            Py_DECREF(_array_differential);
            Py_DECREF(_array_algebraic);
            translate_python_exception();
        }

        const Real retval = PyFloat_AsDouble(_result);

        Py_DECREF(_result);
        Py_DECREF(_current_time);
        Py_DECREF(_array_differential);
        Py_DECREF(_array_algebraic);
        return retval;
    }

private:
    PyObject* callable_;
};

#endif /* PYTHONFUNCTION_HPP */
