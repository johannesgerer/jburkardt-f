#!/bin/bash
#
#  Compile with GFORTRAN.
#
gfortran -fopenmp compute_pi.f90
#
#  Compile with IFORT.
#
#ifort -openmp -parallel -fpp compute_pi.f90
#
mv a.out compute_pi
#
#  Request 1 thread.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./compute_pi > compute_pi_local_output.txt
#
#  Request 2 threads.
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./compute_pi >> compute_pi_local_output.txt
#
#  Request 4 threads.
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./compute_pi >> compute_pi_local_output.txt
#
#  Discard the executable.
#
rm compute_pi
#
echo "Program output written to compute_pi_local_output.txt"