#!/bin/bash
#
gfortran -c -g f90_intrinsics.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling f90_intrinsics.f90"
  exit
fi
rm compiler.txt
#
gfortran f90_intrinsics.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading f90_intrinsics.o"
  exit
fi
rm f90_intrinsics.o
#
mv a.out f90_intrinsics
./f90_intrinsics > f90_intrinsics_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running f90_intrinsics"
  exit
fi
rm f90_intrinsics
#
echo "Program output written to f90_intrinsics_output.txt."
