#!/bin/bash
#
#  Compile with GFORTRAN.
#
gfortran -fopenmp quad2d_openmp.f90
#
mv a.out quad2d
#
#  Request 1 thread.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./quad2d > quad2d_local_gfortran_output.txt
#
#  Request 2 threads.
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./quad2d >> quad2d_local_gfortran_output.txt
#
#  Request 4 threads.
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./quad2d >> quad2d_local_gfortran_output.txt
#
#  Discard the executable.
#
rm quad2d
#
echo "Program output written to quad2d_local_gfortran_output.txt"
