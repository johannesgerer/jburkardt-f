#!/bin/bash
#
gfortran -c -g frieze.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling frieze.f90"
  exit
fi
rm compiler.txt
#
gfortran frieze.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading frieze.o"
  exit
fi
#
rm frieze.o
mv a.out ~/bin/$ARCH/frieze
#
echo "Program installed as ~/bin/$ARCH/frieze"
