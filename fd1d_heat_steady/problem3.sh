#!/bin/bash
#
gfortran-c -g problem3.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling problem3.f90"
  exit
fi
rm compiler.txt
#
gfortran problem3.o -L$HOME/lib/$ARCH -lfd1d_heat_steady
if [ $? -ne 0 ]; then
  echo "Errors linking and loading problem3.o"
  exit
fi
rm problem3.o
#
mv a.out problem3
./problem3 > problem3_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running problem3"
  exit
fi
rm problem3
#
echo "Test program output written to problem3_output.txt."
