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
gfortran -fopenmp hello_openmp.f90
#
mv a.out hello
#
#  Run with 1, 2, 4 and 8 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./hello > hello_ithaca_gfortran_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./hello >> hello_ithaca_gfortran_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./hello >> hello_ithaca_gfortran_output.txt
#
echo "Run with 8 threads."
export OMP_NUM_THREADS=8
./hello >> hello_ithaca_gfortran_output.txt
#
#  Clean up.
#
rm hello
#
echo "Program output written to hello_ithaca_gfortran_output.txt."
