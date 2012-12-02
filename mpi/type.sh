#!/bin/bash
#
gfortran -c type.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling type.f90"
  exit
fi
rm compiler.txt
#
gfortran type.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading type.o"
  exit
fi
rm type.o
#
mv a.out type
mpirun -v -np 4 ./type > type_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running type"
  exit
fi
rm type
#
echo "Program output written to type_output.txt"
