convolve: convolve.o
	g++ -o convolve convolve.o

convolve.o: convolve.cpp
	g++ -c convolve.cpp