from unittest import TestCase, main
import model

class ModelTest(TestCase):
    def setUp(self):
        self.m = model.Model()
        phosphorylation_state = model.StateType(
            'phosphorylated', [
                'unphosphorylated',
                'phosphorylated'])

        self.S1 = self.m.new_species_type(
            name='S1',
            states={'p': phosphorylation_state})
        self.S2 = self.m.new_species_type(name='S2')
        self.P = self.m.new_species_type(name='P')

    def test_basic(self):
        # reversible reaction
        rules = self.m.network_rules
        rules.add([self.S1, self.S2], [self.P], k=.1)
        rules.add([self.P], [self.S1, self.S2], k=.3)

        S1_p = rules.get_species(self.S1, {'p': 'phosphorylated'})

        rr = rules.query(self.S1, self.S2)
        self.assertEqual(rr['k'], .1)

        rr = rules.query(self.P)
        self.assertEqual(rr['k'], .3)
        rules.add([S1_p, self.S2], [self.P], k=.2)

        rr = rules.query(S1_p, self.S2)
        self.assertEqual(rr['k'], .2)

    def test_all_species(self):
        # reversible reaction
        rules = self.m.network_rules
        rules.add([self.S1, self.S2], [self.P], k=.1)
        rules.add([self.P], [self.S1, self.S2], k=.3)

        all_species = set(rules.all_species())
        self.assertEqual(len(all_species), 5)

    def test_generate_all(self):
        rules = self.m.network_rules
        rules.add([self.S1, self.S2], [self.P], k=.1)
        rules.add([self.P], [self.S1, self.S2], k=.3)

        all_rules = set(rules.generate_all())
        self.assertEqual(len(all_rules), 4)
        

if __name__ == '__main__':
    main()
