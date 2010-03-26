class StateType(object):
    def __init__(self, id, states=[]):
        self.__id = id
        self.__states = ['*'] + list(states)

    def add(self, name):
        self.__states.append(name)

    def get_name(self, id):
        return self.__states[id]

    def get_value(self, name):
        for i, n in enumerate(self.__states):
            if n == name:
                return i

    def itervalues(self):
        return xrange(1, len(self.__states))

    def iternames(self):
        return self.__states[1:]

    def __str__(self):
        return str(self.__states)


class SpeciesType(object):
    def __init__(self, id=None):
        self.id = id
        self.__attrs = {}
        self.__state_types = {}

    def add_state_type(self, name, state_type):
        self.__state_types[name] = state_type

    def get_state_type(self, name):
        return self.__state_types[name]

    @property
    def state_types(self):
        return self.__state_types.iteritems()

    @property
    def state_type_count(self):
        return len(self.__state_types)

    def __setitem__(self, k, v):
        self.__attrs[k] = v

    def __getitem__(self, k):
        return self.__attrs[k]

    @property
    def attributes(self):
        return self.__attrs.iteritems()

    def __cmp__(self, rhs):
        return cmp(id(self), id(rhs))

    def __hash__(self):
        return id(self)

    def __str__(self):
        retval = 'SpeciesType('
        retval += 'id=%s, ' % self.id
        retval += 'state_types={'
        for i, (name, state_type) in enumerate(self.__state_types.iteritems()):
            if i > 0:
                retval += ', '
            retval += '\'%s\': %s' % (name, state_type)
        retval += '}, '
        retval += 'attrs=%s' % self.__attrs
        retval += ')'
        return retval


class Species(object):
    def __init__(self, species_type, id=None):
        self.__species_type = species_type
        self.id = id
        states = {}
        for name, _ in species_type.state_types:
            states[name] = 0
        self.__states = states
        self.__attrs = {}

    @property
    def species_type(self):
        return self.__species_type

    def set_state(self, name, value):
        if isinstance(value, str):
            value = self.species_type.get_state_type(name).get_value(value)
        self.__states[name] = value

    def get_state(self, name):
        return self.__states[name]

    def __setitem__(self, k, v):
        self.__attrs[k] = v

    def __getitem__(self, k):
        return self.__attrs[k]

    @property
    def attributes(self):
        return self.__attrs.iteritems()

    def __cmp__(self, rhs):
        a = cmp(self.__species_type, rhs.__species_type)
        if a != 0:
            return a
        return cmp(self.__states, rhs.__states)

    def __eq__(self, rhs):
        return cmp(self, rhs) == 0

    def __hash__(self):
        retval = hash(self.__species_type)
        for v in self.__states.iteritems():
            retval ^= hash(v)
        return retval

    def __str__(self):
        retval = 'Species('
        retval += 'id=%s, ' % self.id
        retval += 'species_type=%s, ' % self.__species_type
        retval += 'states={'
        for i, (name, state_value) in enumerate(self.__states.iteritems()):
            state_type = self.__species_type.get_state_type(name)
            if i > 0:
                retval += ', '
            retval += '\'%s\': \'%s\'' % (name, state_type.get_name(state_value))
        retval += '}, '
        retval += 'attrs=%s' % self.__attrs
        retval += ')'
        return retval


class ReactionRule(object):
    def __init__(self, reactants, products):
        self.__reactants = tuple(sorted(reactants))
        self.__products = tuple(sorted(products))
        self.__attrs = {}

    @property
    def reactants(self):
        return self.__reactants

    @property
    def products(self):
        return self.__products

    def __setitem__(self, k, v):
        self.__attrs[k] = v

    def __getitem__(self, k):
        return self.__attrs.get(k, None)

    @property
    def attributes(self):
        return self.__attrs.iteritems()

    def __cmp__(self, rhs):
        a = cmp(self.__reactants, rhs.__reactants)
        if a != 0:
            return a
        return cmp(self.__products, self.__products)

    def __str__(self):
        retval = 'ReactionRule('
        retval += 'reactants=('
        for i, species in enumerate(self.__reactants):
            if i > 0:
                retval += ', '
            retval += str(species)
        retval += '), '
        retval += 'products=('
        for i, species in enumerate(self.__products):
            if i > 0:
                retval += ', '
            retval += str(species)
        retval += '), '
        retval += 'attrs=%r' % self.__attrs
        retval += ')'
        return retval


class NetworkRules(object):
    def __init__(self):
        self.__reaction_rules = {}
        self.__instantiated_species = {}
        self.__serial_to_species_map = {}
        self.__serial = 0

    def get_species(self, species_type, states={}):
        _states = {}
        species = Species(species_type)
        for name, value in states.iteritems():
            species.set_state(name, value)
        serial = self.__instantiated_species.get(species, 0)
        if serial == 0:
            self.__serial += 1
            serial = self.__serial
            species.id = serial
            self.__instantiated_species[species] = serial
            self.__serial_to_species_map[serial] = species
        else:
            species = self.__serial_to_species_map[serial]
        return species

    def __iter__(self):
        return self.__reaction_rules.itervalues()

    def __ensure_instantiated(self, species_list):
        retval = []
        for species in species_list:
            if isinstance(species, SpeciesType):
                species = self.get_species(species)
            retval.append(species)
        return retval

    def add(self, reactants, products, **attrs):
        reaction_rule = ReactionRule(
            self.__ensure_instantiated(reactants),
            self.__ensure_instantiated(products))
        self.__reaction_rules[reaction_rule.reactants] = reaction_rule
        for k, v in attrs.iteritems():
            reaction_rule[k] = v
        return reaction_rule

    def query(self, s1, s2=None):
        if s2 is None:
            reactants = self.__ensure_instantiated([s1])
            return self.__reaction_rules.get(tuple(reactants), None)
        else:
            reactants = self.__ensure_instantiated([s1, s2])
            return self.__reaction_rules.get(tuple(reactants), None)

    @classmethod
    def __enumerate_all_states(self, states, state_type_list, skip_variant=False):
        if not state_type_list:
            yield states
            return
        name, state_type = state_type_list.pop()
        if not skip_variant:
            states[name] = 0
            for i in self.__enumerate_all_states(states, state_type_list):
                yield i
        for state_value in state_type.itervalues():
            states[name] = state_value
            for i in self.__enumerate_all_states(states, state_type_list):
                yield i

    def all_species(self):
        all_species_types = set()
        for reaction_rule in self.__reaction_rules.itervalues():
            for reactant in reaction_rule.reactants:
                all_species_types.add(reactant.species_type)
            for product in reaction_rule.products:
                all_species_types.add(product.species_type)
        for species_type in all_species_types:
            for states in self.__enumerate_all_states({}, list(species_type.state_types)):
                yield self.get_species(species_type, states)

    def __generate_combinations(self, species_list, species_types):
        if not species_types:
            yield species_list
            return

        species_type = species_types.pop()
        idx = len(species_list)
        species_list.append(None)

        for states in self.__enumerate_all_states({}, list(species_type.state_types), True):
            species_list[idx] = self.get_species(species_type, states)
            yield species_list

    def __enumerate_invariant_combinations(self, species_list):
        invariant_species = []
        variant_species_types = []

        for species in species_list:
            if species.species_type.state_type_count == 0:
                invariant_species.append(species)
            else:
                for name, state_type in species.species_type.state_types:
                    if species.get_state(name) == 0:
                        variant_species_types.append(species.species_type)
                    else:
                        invariant_species.append(species)

        for i in self.__generate_combinations([], variant_species_types):
            yield invariant_species + i

    def generate_all(self):
        for reaction_rule in self.__reaction_rules.itervalues():
            for reactants in self.__enumerate_invariant_combinations(reaction_rule.reactants):
                for products in self.__enumerate_invariant_combinations(reaction_rule.products):
                    rr = ReactionRule(reactants, products)
                    for k, v in reaction_rule.attributes:
                        rr[k] = v
                    yield rr


class Model(object):
    def __init__(self):
        self.__species_types = []
        self.__network_rules = NetworkRules()
        self.__serial = 0

    def new_species_type(self, states={}, **attrs):
        self.__serial += 1
        species_type = SpeciesType(self.__serial)
        for name, state_type in states.iteritems():
            assert(isinstance(state_type, StateType))
            species_type.add_state_type(name, state_type)
        for k, v in attrs.iteritems():
            species_type[k] = v
        self.__species_types.append(species_type)
        return species_type

    @property
    def network_rules(self):
        return self.__network_rules
