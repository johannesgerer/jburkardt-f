#!/bin/bash
#
cp ~/include/mpi_stubs_f90.h mpif.h
#
gfortran -c -g quadrature.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors occurred while compiling quadrature.f90"
  exit
fi
rm compiler.txt
rm mpif.h
#
gfortran quadrature.o -L$HOME/lib/$ARCH -lmpi_stubs
if [ $? -ne 0 ]; then
  echo "Errors occurred while linking and loading quadrature.o"
  exit
fi
rm quadrature.o
#
mv a.out quadrature
./quadrature > quadrature_output.txt
if [ $? -ne 0 ]; then
  echo "Errors occurred while running quadrature"
  exit
fi
rm quadrature
#
echo "Program output written to quadrature_output.txt"
