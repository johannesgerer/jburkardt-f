#!/bin/bash
#
#  Compile with IFORT.
#
ifort -openmp -parallel -fpp quad_openmp.f90
#
mv a.out quad
#
#  Request 1 thread.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./quad > quad_local_ifort_output.txt
#
#  Request 2 threads.
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./quad >> quad_local_ifort_output.txt
#
#  Request 4 threads.
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./quad >> quad_local_ifort_output.txt
#
#  Discard the executable.
#
rm quad
#
echo "Program output written to quad_local_ifort_output.txt"
