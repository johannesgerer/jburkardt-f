#!/bin/bash
#
#  Compile with IFORT.
#
ifort -openmp -parallel -fpp fft_openmp.f90
#
mv a.out fft
#
#  Request 1 thread.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./fft > fft_local_ifort_output.txt
#
#  Request 2 threads.
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./fft >> fft_local_ifort_output.txt
#
#  Request 4 threads.
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./fft >> fft_local_ifort_output.txt
#
#  Discard the executable.
#
rm fft
#
echo "Program output written to fft_local_ifort_output.txt"
