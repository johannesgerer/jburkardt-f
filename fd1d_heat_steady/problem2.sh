#!/bin/bash
#
gfortran -c -g problem2.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling problem2.f90"
  exit
fi
rm compiler.txt
#
gfortran problem2.o -L$HOME/lib/$ARCH -lfd1d_heat_steady
if [ $? -ne 0 ]; then
  echo "Errors linking and loading problem2.o"
  exit
fi
rm problem2.o
#
mv a.out problem2
./problem2 > problem2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running problem2"
  exit
fi
rm problem2
#
echo "Test program output written to problem2_output.txt."
