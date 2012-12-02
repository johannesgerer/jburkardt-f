#!/bin/bash
#
gfortran -c intervals_mpi.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling intervals_mpi.f90"
  exit
fi
rm compiler.txt
#
gfortran intervals_mpi.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading intervals_mpi.o"
  exit
fi
rm intervals_mpi.o
#
mv a.out intervals
mpirun -v -np 4 ./intervals > intervals_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running intervals"
  exit
fi
rm intervals
#
echo "Program output written to intervals_output.txt"
