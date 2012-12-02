#!/bin/bash
#
gfortran -c -g gfortran_intrinsics.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling gfortran_intrinsics.f90"
  exit
fi
rm compiler.txt
#
gfortran gfortran_intrinsics.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading gfortran_intrinsics.o"
  exit
fi
rm gfortran_intrinsics.o
#
mv a.out gfortran_intrinsics
./gfortran_intrinsics > gfortran_intrinsics_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running gfortran_intrinsics"
# exit
fi
rm gfortran_intrinsics
#
echo "Program output written to gfortran_intrinsics_output.txt."
