LIBS = -L/usr/lib -lgsl -lgslcblas -lm

all: main

main: main.cpp DAESolver.cpp
	g++ -o main -O2 $(LIBS) *.cpp

clean:
	rm main

