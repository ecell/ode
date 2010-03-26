#ifndef UTILS_HPP
#define UTILS_HPP

#include <stdexcept>
#include <string>
#include <Python.h>

class PythonException: public std::exception
{
public:
    PythonException(PyObject* errType, PyObject* errValue,
        PyObject* traceback)
      :
      errType_(errType), errValue_(errValue), traceback_(traceback_)
    {
        Py_INCREF(errType_);
        Py_INCREF(errValue_);
        Py_INCREF(traceback_);
    }

    virtual ~PythonException() throw()
    {
        Py_DECREF(traceback_);
        Py_DECREF(errValue_);
        Py_DECREF(errType_);
    }

    void restore()
    {
        PyErr_Restore(errType_, errValue_, traceback_);
    }

    virtual const char* what() const throw()
    {
        if (whatmsg_.empty())
        {
            char *msg = 0;
            Py_ssize_t msglen = 0;

            PyObject* msgobj = PyObject_Str(errValue_);
            if (msgobj)
            {
                PyString_AsStringAndSize(msgobj, &msg, &msglen);
                whatmsg_.assign(Py_TYPE(errValue_)->tp_name);
                if (msglen)
                {
                    whatmsg_.insert(whatmsg_.size(), ": ");
                    whatmsg_.insert(whatmsg_.size(), msg, msglen);
                }
                Py_DECREF(msgobj);
            }
            else
            {
                whatmsg_.assign(Py_TYPE(errValue_)->tp_name);
            }
        }
        return whatmsg_.c_str();
    }

private:
    PyObject* errType_;
    PyObject* errValue_;
    PyObject* traceback_;
    mutable std::string whatmsg_;
};

inline void translate_python_exception()
{
    PyObject* errType = 0, *errValue = 0, *traceback = 0;
    PyErr_Fetch(&errType, &errValue, &traceback);
    PyErr_Clear();
    throw PythonException(errType, errValue, traceback);
}

#endif /* UTILS_HPP */
