#!/bin/bash
#
cp ~/include/mpi_stubs_f90.h mpif.h
#
gfortran -c -g buffon_laplace.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors occurred while compiling buffon_laplace.f90"
  exit
fi
rm compiler.txt
rm mpif.h
#
gfortran buffon_laplace.o -L$HOME/lib/$ARCH -lmpi_stubs
if [ $? -ne 0 ]; then
  echo "Errors occurred while linking and loading buffon_laplace.o"
  exit
fi
rm buffon_laplace.o
#
mv a.out buffon_laplace
./buffon_laplace > buffon_laplace_output.txt
if [ $? -ne 0 ]; then
  echo "Errors occurred while running buffon_laplace"
  exit
fi
rm buffon_laplace
#
echo "Program output written to buffon_laplace_output.txt"
