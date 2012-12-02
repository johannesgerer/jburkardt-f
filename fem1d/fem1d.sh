#!/bin/bash
#
gfortran -c -g fem1d.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling fem1d.f90"
  exit
fi
rm compiler.txt
#
gfortran fem1d.o
if [ $? -ne 0 ]; then
  echo "Errors while loading fem1d.o"
  exit
fi
rm fem1d.o
#
mv a.out ~/bin/$ARCH/fem1d
#
echo "Program installed as ~/bin/$ARCH/fem1d"
