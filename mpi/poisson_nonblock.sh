#!/bin/bash
#
gfortran -c poisson_nonblock_mpi.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling poisson_nonblock_mpi.f90"
  exit
fi
rm compiler.txt
#
gfortran poisson_nonblock_mpi.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading poisson_nonblock_mpi.o"
  exit
fi
rm poisson_nonblock_mpi.o
#
mv a.out poisson_nonblock
mpirun -v -np 4 ./poisson_nonblock > poisson_nonblock_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running poisson_nonblock"
  exit
fi
rm poisson_nonblock
#
echo "Program output written to poisson_nonblock_output.txt"
