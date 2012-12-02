#!/bin/bash
#
gfortran -c mxm.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mxm.f90"
  exit
fi
rm compiler.txt
#
gfortran mxm.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mxm.o"
  exit
fi
rm mxm.o
#
mv a.out ~/bin/$ARCH/mxm
#
echo "Executable installed as ~/bin/$ARCH/mxm"
