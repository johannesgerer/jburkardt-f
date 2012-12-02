#!/bin/bash
#
gfortran -c -g problem1.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling problem1.f90"
  exit
fi
rm compiler.txt
#
gfortran problem1.o -L$HOME/lib/$ARCH -lfd1d_heat_steady
if [ $? -ne 0 ]; then
  echo "Errors linking and loading problem1.o"
  exit
fi
rm problem1.o
#
mv a.out problem1
./problem1 > problem1_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running problem1"
  exit
fi
rm problem1
#
echo "Test program output written to problem1_output.txt."
