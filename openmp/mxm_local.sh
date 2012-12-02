#!/bin/bash
#
#  Compile with GFORTRAN.
#
gfortran -fopenmp mxm.f90
#
#  Compile with IFORT.
#
#ifort -openmp -parallel -fpp mxm.f90
#
mv a.out mxm
#
#  Request 1 thread.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./mxm > mxm_local_output.txt
#
#  Request 2 threads.
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./mxm >> mxm_local_output.txt
#
#  Request 4 threads.
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./mxm >> mxm_local_output.txt
#
#  Discard the executable.
#
rm mxm
#
echo "Program output written to mxm_local_output.txt"