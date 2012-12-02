#!/bin/bash
#
gfortran -c -g fem1d_pmethod.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling fem1d_pmethod.f90"
  exit
fi
rm compiler.txt
#
gfortran fem1d_pmethod.o
if [ $? -ne 0 ]; then
  echo "Errors while loading fem1d_pmethod.o"
  exit
fi
rm fem1d_pmethod.o
#
mv a.out ~/bin/$ARCH/fem1d_pmethod
#
echo "Program installed as ~/bin/$ARCH/fem1d_pmethod"
