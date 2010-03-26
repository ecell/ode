import numpy
from process.process import FunctionMaker

class Simulator:
    def __init__(self):
        self.the_time = 0.0
        self.dimension = 0
        self.dimension_algebraic = 0
        self.outputs_series = []

    def make_functions(self, model, volume):
        '''Make functions from model
        '''
        function_maker = FunctionMaker()
        return function_maker.make_functions(model, volume)

    def initialize(self, a_solver, functions, variables,
            variables_algebraic = []):
        '''Initialize solver

        1. Set variables.
        2. Initialize global time.
        3. Clear outputs.
        4. Set functions.
        '''
        if len(functions) != len(variables) + len(variables_algebraic):
            print 'Number of function and number \
                of variable have to be equal.'
            return -1

        self.solver = a_solver

        self.dimension = len(variables)
        variable_array = numpy.array(variables)
        if variables_algebraic != []:
            # using DAESolver
            self.dimension_algebraic = len(variables_algebraic)
            self.dimension += self.dimension_algebraic
            variable_algebraic_array = numpy.array(variables_algebraic)
            self.solver.initialize(variable_array, variable_algebraic_array)
        else:
            self.dimension_algebraic = 0
            self.solver.initialize(variable_array)

        self.the_time = self.solver.get_time()
        self.outputs_series = []

        # Register Functions
        for function in functions:
            self.__register_function(function)

        return 0

    def __register_function(self, function):
        '''Set a first order ordinal diferential equation.'''
        self.solver.register_function(function)

    def __log_data(self):
        '''Log the local time and the variables.'''
        outputs = numpy.array([self.solver.get_time()])
        variable_array = numpy.empty(self.dimension)

        if self.dimension_algebraic > 0:
            variable_differential_array = numpy.empty(
                self.dimension - self.dimension_algebraic)
            variable_algebraic_array = numpy.empty(self.dimension_algebraic)
            self.solver.get_variable_differential_array(
                variable_differential_array)
            self.solver.get_variable_algebraic_array(
                variable_algebraic_array)
            variable_array = numpy.append(variable_differential_array,
                variable_algebraic_array)
        else:
            self.solver.get_variable_array(variable_array)

        outputs = numpy.append(outputs, variable_array)
        self.outputs_series.append(outputs)

    def get_logged_data(self):
        '''Get the logged data series.'''
        return numpy.array(self.outputs_series)

    def get_time(self):
        '''Get the current time.'''
        return self.solver.get_time()

    def get_variable(self, index):
        '''Get a variable value.'''
        return self.solver.get_variable(index)

    def _step(self):
        self.solver.integrate(self.the_time)
        state_event = self.solver.step()
        self.the_time = self.solver.reschedule()
        self.__log_data()
        return state_event

    def step(self, n = 1):
        '''Solve n step.'''
        try:
            for i in range(n):
                if self._step() != 0:
                    break
        except:
            pass

    def run(self, an_end_time = -1.0):
        '''Run the solver'''
        try:
            while self.the_time < an_end_time:
                if self._step() != 0:
                    break
        except:
            pass

    def register_status_event(self, status_event_list):
        '''Register status events.

        status_event_list = [status event 1, ..., status event N]
            status event is a dictionary register following attribute.
            'id': index of variable array.
            'thres': if variable value is reached this value,
                     the event is araised.
            'code': status code for this status event.
            'flag': 0 indicate the variables_differential,
                    1 indicate the variables_algebraic.
        '''
        for status_event in status_event_list:
            self.solver.register_status_event(status_event)

