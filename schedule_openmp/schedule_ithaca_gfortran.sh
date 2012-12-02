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
gfortran -fopenmp schedule_openmp.f90
#
mv a.out schedule
#
#  Run with 2 threads.
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./schedule > schedule_ithaca_gfortran_output.txt
#  Clean up.
#
rm schedule
#
echo "Program output written to schedule_ithaca_gfortran_output.txt."
