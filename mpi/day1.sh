#!/bin/bash
#
gfortran -c day1.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling day1.f90"
  exit
fi
rm compiler.txt
#
gfortran day1.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading day1.o"
  exit
fi
rm day1.o
#
mv a.out day1
mpirun -v -np 4 ./day1 > day1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running day1"
  exit
fi
rm day1
#
echo "Program output written to day1_output.txt"
