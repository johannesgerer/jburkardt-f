#!/bin/bash
#
gfortran -c mxv.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mxv.f90"
  exit
fi
rm compiler.txt
#
gfortran mxv.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mxv.o"
  exit
fi
rm mxv.o
#
mv a.out ~/bin/$ARCH/mxv
#
echo "Executable installed as ~/bin/$ARCH/mxv"
