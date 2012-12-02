#!/bin/bash
#
#  Compile with GFORTRAN.
#
gfortran -fopenmp hello_openmp.f90
#
mv a.out hello
#
#  Request 1 thread.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./hello > hello_local_gfortran_output.txt
#
#  Request 2 threads.
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./hello >> hello_local_gfortran_output.txt
#
#  Request 4 threads.
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./hello >> hello_local_gfortran_output.txt
#
#  Request 8 threads.
#
echo "Run with 8 threads."
export OMP_NUM_THREADS=8
./hello >> hello_local_gfortran_output.txt
#
#  Discard the executable.
#
rm hello
#
echo "Program output written to hello_local_gfortran_output.txt"
