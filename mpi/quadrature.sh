#!/bin/bash
#
gfortran -c quadrature.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quadrature.f90"
  exit
fi
rm compiler.txt
#
gfortran quadrature.o -lmpi
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quadrature.o"
  exit
fi
rm quadrature.o
#
mv a.out quadrature
mpirun -v -np 4 ./quadrature > quadrature_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quadrature"
  exit
fi
rm quadrature
#
echo "Program output written to quadrature_output.txt"
