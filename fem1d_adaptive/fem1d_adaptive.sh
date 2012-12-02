#!/bin/bash
#
gfortran -c -g fem1d_adaptive.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling fem1d_adaptive.f90"
  exit
fi
rm compiler.txt
#
gfortran fem1d_adaptive.o
if [ $? -ne 0 ]; then
  echo "Errors while loading fem1d_adaptive.o"
  exit
fi
rm fem1d_adaptive.o
#
mv a.out ~/bin/$ARCH/fem1d_adaptive
#
echo "Program installed as ~/bin/$ARCH/fem1d_adaptive"
