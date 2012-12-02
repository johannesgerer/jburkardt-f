#!/bin/bash
#
cp ~/include/mpi_stubs_f90.h mpif.h
#
gfortran -c -g hello.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors occurred while compiling hello.f90"
  exit
fi
rm compiler.txt
rm mpif.h
#
gfortran hello.o -L$HOME/lib/$ARCH -lmpi_stubs
if [ $? -ne 0 ]; then
  echo "Errors occurred while linking and loading hello.o"
  exit
fi
rm hello.o
#
mv a.out hello
./hello > hello_output.txt
if [ $? -ne 0 ]; then
  echo "Errors occurred while running hello"
  exit
fi
rm hello
#
echo "Program output written to hello_output.txt"
