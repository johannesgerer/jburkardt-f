#!/bin/bash
#
#PBS -lwalltime=00:05:00
#PBS -lnodes=1:ppn=8
#PBS -W group_list=ithaca
#PBS -q ithaca_q
#PBS -A admn0000
#PBS -j oe
#
#  Choose the MPI distribution.
#
#export MPIDISTRO=mvapich2-1.4
export MPIDISTRO=openmpi-1.3.3

export LD_LIBRARY_PATH=/apps/packages/$MPIDISTRO/lib:$LD_LIBRARY_PATH

export PATH=/apps/packages/$MPIDISTRO/bin:$PATH

cd $PBS_O_WORKDIR
#
#  Compile
#
mpif90 hello_mpi.f90
mv a.out hello

mpirun -np 8 ./hello > hello_sysx_output.txt
#
rm hello
#
echo "Program output written to hello_sysx_output.txt."
