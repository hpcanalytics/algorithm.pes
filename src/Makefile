#CC = xlc_r
CPP = g++
CFLAGS = -O2 -Wno-deprecated  #-qsmp #-O2 #-qarch=auto -qtune=auto  

all: pshuffle
	
pshuffle: Parallel_Shuffle.cpp Parallel_Shuffle.hpp Timer.hpp Serial_Multinomial.hpp utility.hpp utility2.hpp TSet.hpp Random.hpp
	module load mpi/mvapich2/gcc/4.7.2/1.9-psm;	mpicxx -o pshuffle Parallel_Shuffle.cpp

clean:
	rm -f pshuffle
	

