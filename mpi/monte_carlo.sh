#!/bin/bash
#
gfortran -c monte_carlo.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling monte_carlo.f90"
  exit
fi
rm compiler.txt
#
gfortran monte_carlo.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading monte_carlo.o"
  exit
fi
rm monte_carlo.o
#
mv a.out monte_carlo
mpirun -v -np 4 ./monte_carlo > monte_carlo_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running monte_carlo"
  exit
fi
rm monte_carlo
#
echo "Program output written to monte_carlo_output.txt"
