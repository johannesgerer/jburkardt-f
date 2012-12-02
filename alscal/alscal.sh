#!/bin/bash
#
gfortran -c alscal.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling alscal.f90"
  exit
fi
rm compiler.txt
#
gfortran alscal.o
if [ $? -ne 0 ]; then
  echo "Errors loading alscal.o"
  exit
fi
#
rm alscal.o
mv a.out ~/bin/$ARCH/alscal
#
echo "Executable installed as ~/bin/$ARCH/alscal"
