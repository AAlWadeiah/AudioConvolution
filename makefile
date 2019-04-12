convolve.out: convolve.o
	g++ -o convolve.out convolve.o

convolve.o: convolve.cpp
	g++ -c convolve.cpp