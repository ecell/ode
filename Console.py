#!/usr/bin/env python

import sys
import os
import code
from Simulator import Simulator

class Prompt:
    def __init__(self, a_simulator):
        self.simulator = a_simulator

    def __str__(self):
        return '<t=%g>>> ' %\
            self.simulator.get_time()

class Console:
    def __init__(self, a_simulator, banner):
        self.simulator = a_simulator
        self.banner = banner

    def __create_script_context(self):
        # the_simulator == self in the script
        a_context = {'the_simulator': self.simulator, 'self': self.simulator}

        a_key_list = dir(self.simulator)
        a_dict = {}
        for a_key in a_key_list:
            if not a_key.startswith('__'):
                a_dict[a_key] = getattr(self.simulator, a_key)

        a_context.update(a_dict)
       
        return a_context

    # Session methods
    def load_script(self, ecs):
        (ecs_dir, ecs_file) = os.path.split(ecs)
        if ecs_dir != '':
            os.chdir(ecs_dir)
        a_context = self.__create_script_context()
        execfile(ecs_file, a_context)

    def interact(self):
        a_context = self.__create_script_context()
        try:
            import readline # to provide convenient commandline editing :)
        except:
            pass
        a_context['__prompt__'] = Prompt(self.simulator)
        an_interpreter = code.InteractiveConsole(a_context)
        an_interpreter.runsource(
            'import sys; sys.ps1 = __prompt__; del sys, __prompt__')
        an_interpreter.interact(self.banner)

def main():
    an_ess_file = None

    # -------------------------------------
    # when ecell3-session [essfile] check file exist
    # and set anESSfile 
    # -------------------------------------

    for i in range(len(sys.argv) - 1):
        i = i + 1
        if not sys.argv[i][0] == "-":
            if os.path.isfile(sys.argv[i]):
                if sys.argv[i - 1] != "-e" and sys.argv[i - 1] != "-f":
                    an_ess_file = sys.argv[i]
            else:
                sys.stderr.write("Error: %s does not exist.\n" % sys.argv[i])
                sys.exit(1) 

    a_simulator = Simulator()
    a_console = Console(a_simulator, '')

    if an_ess_file:
        a_console.load_script(an_ess_file)
    else:
        a_console.interact()

if __name__ == '__main__':
    main()
    sys.exit(0)
