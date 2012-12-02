#!/bin/bash
#
gfortran -c -g sparsekit_prb07.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparsekit_prb07.f90"
  exit
fi
rm compiler.txt
#
gfortran sparsekit_prb07.o -L$HOME/lib/$ARCH -lsparsekit
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparsekit_prb07.o"
  exit
fi
rm sparsekit_prb07.o
#
mv a.out sparsekit_prb07
./sparsekit_prb07 > sparsekit_prb07_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparsekit_prb07"
  exit
fi
rm sparsekit_prb07
#
echo "Program output written to sparsekit_prb07_output.txt"
