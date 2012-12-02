#!/bin/bash
#
#PBS -lwalltime=00:10:00
#PBS -lnodes=1:ppn=8
#PBS -W group_list=ithaca
#PBS -q ithaca_q
#PBS -A admn0000
#PBS -j oe
#
#  Start in the directory from which this job was submitted.
#
cd $PBS_O_WORKDIR
#
#  Compile the program with GFORTRAN.
#
gfortran -fopenmp fft_openmp.f90
#
mv a.out fft
#
#  Run with 1, 2, 4 and 8 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./fft > fft_ithaca_gfortran_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./fft >> fft_ithaca_gfortran_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./fft >> fft_ithaca_gfortran_output.txt
#
echo "Run with 8 threads."
export OMP_NUM_THREADS=8
./fft >> fft_ithaca_gfortran_output.txt
#
#  Clean up.
#
rm fft
#
echo "Program output written to fft_ithaca_gfortran_output.txt."
