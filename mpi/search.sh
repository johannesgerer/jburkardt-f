#!/bin/bash
#
gfortran -c search.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling search.f90"
  exit
fi
rm compiler.txt
#
gfortran search.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading search.o"
  exit
fi
rm search.o
#
mv a.out search
mpirun -v -np 4 ./search > search_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running search"
  exit
fi
rm search
#
echo "Program output written to search_output.txt"
