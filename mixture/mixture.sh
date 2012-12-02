#!/bin/bash
#
gfortran -c mixture.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mixture.f90"
  exit
fi
rm compiler.txt
#
gfortran mixture.o
if [ $? -ne 0 ]; then
  echo "Errors loading mixture.o"
  exit
fi
#
rm mixture.o
mv a.out ~/bin/$ARCH/mixture
#
echo "Executable installed as ~/bin/$ARCH/mixture"
