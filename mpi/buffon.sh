#!/bin/bash
#
gfortran -c buffon_mpi.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling buffon_mpi.f90"
  exit
fi
rm compiler.txt
#
gfortran buffon_mpi.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading buffon_mpi.o"
  exit
fi
rm buffon_mpi.o
#
mv a.out buffon
mpirun -v -np 4 ./buffon > buffon_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running buffon"
  exit
fi
rm buffon
#
echo "Program output written to buffon_output.txt"
