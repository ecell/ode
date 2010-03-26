class MassActionFluxProcess:
    def __init__(self, k_value, volume, reactants):
        self.volume = volume
        self.k_value = k_value
        self.reactants = reactants
        self.N_A = 6.0221367e+23

    def __call__(self, variable_array, time):
        velocity = self.k_value * self.volume * self.N_A
        for r in self.reactants:
            coefficient = r['coef']
            value = variable_array[r['id']]
            while coefficient > 0:
                velocity *= value / (self.volume * self.N_A)
                coefficient -= 1
        return velocity

class Function:
    def __init__(self):
        self.process_list = []

    def add_process(self, coefficient, process):
        self.process_list.append(
            {'coef': coefficient, 'func': process})

    def __process(self, variable_array, time):
        for process in self.process_list:
            yield process['coef'] * process['func'](
                variable_array, time)

    def __call__(self, variable_array, time):
        return sum(self.__process(variable_array, time)) + 0.0

class FunctionMaker:
    def __create_rule_list(self, network_rules, dimension):
        '''Create rule list from network rules of model.'''
        rule_list = []
        for r in network_rules.generate_all():
            rule = {}

            rule['k'] = r['k']

            rule['reactants'] = []
            for reactant in r.reactants:
                variable_id = reactant.id - 1
                rule['reactants'].append(
                    {'id': variable_id, 'coef': 1})

            rule['products'] = []
            for product in r.products:
                variable_id = product.id - 1
                rule['products'].append(
                    {'id': variable_id, 'coef': 1})

            rule_list.append(rule)

        return rule_list

    def make_functions(self, model, volume):
        '''Make functions from model
        '''
        dimension = len(set(model.network_rules.all_species()))
        rule_list = self.__create_rule_list(
            model.network_rules, dimension)
        
        #for rule in rule_list:
        #    print rule

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

