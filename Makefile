subsystem:
	cd solver && $(MAKE)
	cd solver/cpp_lib/ODE && $(MAKE)
	cd solver/cpp_lib/DAE && $(MAKE)

clean:
	cd solver && $(MAKE) clean
	cd solver/cpp_lib/ODE && $(MAKE) clean
	cd solver/cpp_lib/DAE && $(MAKE) clean

