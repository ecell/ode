from unittest import TestCase, main
import numpy
from process.process import FunctionMaker
from model.model import Model

class ModelTest(TestCase):
    def setUp(self):
        self.volume = 1e-18
        self.model = Model()
        self.function_maker = FunctionMaker()

    def test_basic(self):
        S1 = self.model.new_species_type(name='S1')
        S2 = self.model.new_species_type(name='S2')
        P = self.model.new_species_type(name='P')

        # S1 + S2 -> P    k=.3
        self.model.network_rules.add([S1, S2], [P], k=.3)
        # P -> S1 + S2    k=.1 (inverse reaction)
        self.model.network_rules.add([P], [S1, S2], k=.1)

        # Make Function list
        functions = self.function_maker.make_functions(self.model, self.volume)

        species = list(self.model.network_rules.all_species())
        self.assertEqual(len(functions), len(species))

if __name__ == '__main__':
    main()
