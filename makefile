convolver.out: convolver.o
	g++ -o convolver.out convolver.o

convolver.o: convolver.cpp
	g++ -c convolver.cpp