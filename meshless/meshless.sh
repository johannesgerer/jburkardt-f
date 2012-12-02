#!/bin/bash
#
gfortran -c -g meshless.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling meshless.f90"
  exit
fi
rm compiler.txt
#
gfortran meshless.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading meshless.o"
  exit
fi
rm meshless.o
#
mv a.out ~/bin/$ARCH/meshless
#
echo "Program installed as ~/bin/$ARCH/meshless"
