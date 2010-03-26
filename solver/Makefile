PYTHON = python

all: Solvers

Solvers: ODESolver.pyx DAESolver.pyx
	$(PYTHON) setup.py build_ext --inplace

clean:
	$(PYTHON) setup.py clean -a
	rm -f ODESolver.cpp
	rm -f ODESolver.so
	rm -f DAESolver.cpp
	rm -f DAESolver.so

.PHOHY: all clean
