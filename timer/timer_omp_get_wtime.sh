#!/bin/bash
#
gfortran -fopenmp timer_omp_get_wtime.f90
mv a.out timer_omp_get_wtime
#
#  Run the program with 1 thread.
#
export OMP_NUM_THREADS=1
#
./timer_omp_get_wtime > timer_omp_get_wtime_output.txt
rm timer_omp_get_wtime
#
echo "Program output written to timer_omp_get_wtime_output.txt"
