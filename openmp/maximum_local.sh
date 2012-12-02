#!/bin/bash
#
#  Compile with GFORTRAN.
#
gfortran -fopenmp maximum.f90
#
#  Compile with IFORT.
#
#ifort -openmp -parallel -fpp maximum.f90
#
mv a.out maximum
#
#  Request 1 thread.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./maximum > maximum_local_output.txt
#
#  Request 2 threads.
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./maximum >> maximum_local_output.txt
#
#  Request 4 threads.
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./maximum >> maximum_local_output.txt
#
#  Discard the executable.
#
rm maximum
#
echo "Program output written to maximum_local_output.txt"