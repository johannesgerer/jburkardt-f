#!/bin/bash
#
#  Compile with GFORTRAN.
#
gfortran -fopenmp helmholtz.f90
#
#  Compile with IFORT.
#
#ifort -openmp -parallel -fpp helmholtz.f90
#
mv a.out helmholtz
#
#  Request 1 thread.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./helmholtz > helmholtz_local_output.txt
#
#  Request 2 threads.
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./helmholtz >> helmholtz_local_output.txt
#
#  Request 4 threads.
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./helmholtz >> helmholtz_local_output.txt
#
#  Discard the executable.
#
rm helmholtz
#
echo "Program output written to helmholtz_local_output.txt"