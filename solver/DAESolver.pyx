__all__ = ["DAESolver"]

cdef extern from "Defs.hpp":
    ctypedef struct StatusEvent:
        long int variable_index
        double threshold
        long int status_code
        bint variable_flag

cdef extern from "Function.hpp":
    ctypedef struct c_Function "Function":
        pass

cdef extern from "PythonFunction.hpp":
    ctypedef struct c_PythonFunction "PythonFunction":
        pass
    c_PythonFunction *new_PythonFunction "new PythonFunction"(object)
    void del_PythonFunction "delete"(c_PythonFunction*) 

cdef extern from "numpy/ndarrayobject.h":
    void import_array() 

cdef extern from "numpy/arrayobject.h":
    ctypedef int intp
    ctypedef extern class numpy.ndarray [object PyArrayObject]:
        cdef char *data
        cdef int nd
        cdef intp *dimensions
        cdef intp *strides
        cdef int flags


cdef extern from "DAESolver.hpp":
    ctypedef struct c_Solver "DAESolver":
        void register_status_event(StatusEvent)
        void register_function(c_Function*)
        void initialize(double*, double*, long int, long int)
        void get_variable_differential_array(double*, long int)
        void get_variable_algebraic_array(double*, long int)
        double get_value_differential(long int)
        double get_value_algebraic(long int)
        double get_current_time()
        long int step()
        void integrate(double)
        double reschedule()
    c_Solver *new_Solver "new DAESolver"()
    void del_Solver "delete"(c_Solver*) 

cdef extern from "CallStep.hpp":
    ctypedef struct StepResult:
        long int state_event
        int err
    StepResult call_step(c_Solver *)

cdef class DAESolver:
    '''Wrapper class of DAESolver'''
    cdef c_Solver *thisptr

    def __cinit__(self):
        self.thisptr = new_Solver()

    def __del__(self):
        del_Solver(self.thisptr)

    def register_function(self, callable):
        '''Set a first order ordinal diferential equation.'''
        cdef c_PythonFunction* pyfun = new_PythonFunction(callable)
        self.thisptr.register_function(<c_Function*>pyfun)

    def initialize(self, ndarray py_variable_differential_array,
            ndarray py_variable_algebraic_array):
        '''Initialize solver
        Set variables.
        '''
        cdef double *variable_differential_array \
            = <double *>py_variable_differential_array.data
        cdef double *variable_algebraic_array \
            = <double *>py_variable_algebraic_array.data

        cdef int dimension_differential \
            = py_variable_differential_array.dimensions[0]
        cdef int dimension_algebraic \
            = py_variable_algebraic_array.dimensions[0]

        self.thisptr.initialize(variable_differential_array,
            variable_algebraic_array,
            dimension_differential + dimension_algebraic,
            dimension_differential)

    def get_variable_differential_array(self, ndarray py_variable_array):
        '''Get the variables for differential.'''
        cdef double *variable_array = <double *>py_variable_array.data
        cdef int dimension = py_variable_array.dimensions[0]
        self.thisptr.get_variable_differential_array(
            variable_array, dimension)

    def get_variable_algebraic_array(self, ndarray py_variable_array):
        '''Get the variables for algebraic.'''
        cdef double *variable_array = <double *>py_variable_array.data
        cdef int dimension = py_variable_array.dimensions[0]
        self.thisptr.get_variable_algebraic_array(variable_array, dimension)

    def get_variable_differential(self, index):
        '''Get the variable value specified by index for differential.'''
        return self.thisptr.get_value_differential(index)

    def get_variable_algebraic(self, index):
        '''Get the variable value specified by index for algebraic.'''
        return self.thisptr.get_value_algebraic(index)

    def get_time(self):
        '''Get the current time.'''
        return self.thisptr.get_current_time()

    def integrate(self, a_time):
        '''Integrate'''
        self.thisptr.integrate(a_time)

    def step(self):
        '''Solve one step.'''
        cdef StepResult step_result
        step_result = call_step(self.thisptr)
        if step_result.err:
            raise Exception
        return step_result.state_event

    def reschedule(self):
        '''Schedule the next time.'''
        return self.thisptr.reschedule()

    def register_status_event(self, py_status_event):
        '''Register a status event'''
        cdef StatusEvent status_event
        status_event.variable_index = py_status_event['id']
        status_event.threshold = py_status_event['thres']
        status_event.status_code = py_status_event['code']
        status_event.variable_flag = py_status_event['flag']
        self.thisptr.register_status_event(status_event)

import_array()
