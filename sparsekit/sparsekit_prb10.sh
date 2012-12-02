#!/bin/bash
#
gfortran -c -g sparsekit_prb10.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparsekit_prb10.f90"
  exit
fi
rm compiler.txt
#
gfortran sparsekit_prb10.o -L$HOME/lib/$ARCH -lsparsekit
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparsekit_prb10.o"
  exit
fi
rm sparsekit_prb10.o
#
mv a.out sparsekit_prb10
./sparsekit_prb10 > sparsekit_prb10_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparsekit_prb10"
  exit
fi
rm sparsekit_prb10
#
echo "Program output written to sparsekit_prb10_output.txt"
