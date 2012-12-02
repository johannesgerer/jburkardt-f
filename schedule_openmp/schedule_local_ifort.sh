#!/bin/bash
#
#  Compile with IFORT.
#
ifort -openmp -parallel -fpp schedule_openmp.f90
#
mv a.out schedule
#
#  Request 2 threads.
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./schedule > schedule_local_ifort_output.txt
#
#  Discard the executable.
#
rm schedule
#
echo "Program output written to schedule_local_ifort_output.txt"
