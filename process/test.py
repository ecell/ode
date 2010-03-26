from unittest import TestCase, main
import numpy
from process import Function, MassActionFluxProcess

def make_functions(rule_list, volume, dimension):
    # Process list
    processes = []
    for rule in rule_list:
        process = MassActionFluxProcess(rule['k'], volume,
            rule['reactants'])
        processes.append(process)

    # Function list
    functions = []
    for i in range(dimension):
        function = Function()
        functions.append(function)

    # Add Processes
    process_id = 0
    for rule in rule_list:
        for reactant in rule['reactants']:
            functions[reactant['id']].add_process(
                -reactant['coef'], processes[process_id])
        for product in rule['products']:
            functions[product['id']].add_process(
                product['coef'], processes[process_id])
        process_id += 1

    return functions


class ModelTest(TestCase):
    def setUp(self):
        self.variable_array = numpy.array([10000., 5000., 10.])
        self.dimension = 3
        self.volume = 1e-18
    def test_basic(self):
        rule_list = [
            {
                'k': 0.1,
                'reactants': [
                    {'id': 0, 'coef': 1},
                    {'id': 1, 'coef': 1}],
                'products': [
                    {'id': 2, 'coef': 1}]
            },
            {
                'k': 0.3,
                'reactants': [
                    {'id': 2, 'coef': 1}],
                'products': [
                    {'id': 0, 'coef': 1},
                    {'id': 1, 'coef': 1}]
            }
        ]

        functions = make_functions(rule_list, self.volume, self.dimension)

        # Calc. velocity
        velocity = []
        for function in functions:
            velocity.append(function(self.variable_array, 0.))
        print velocity

    def test_Error(self):
        rule_list = [
            {
                'k': 0.3,
                'reactants': [
                    {'id': 0, 'coef': 2},
                    {'id': 2, 'coef': 1}],
                'products': [
                    {'id': 1, 'coef': 1}]
            }
        ]

        functions = make_functions(rule_list, self.volume, self.dimension)

        variable_array = numpy.array([10000., 5000.])
        self.assertRaises(IndexError, functions[0],
            variable_array, 0.)

if __name__ == '__main__':
    main()
