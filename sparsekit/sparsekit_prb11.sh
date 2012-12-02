#!/bin/bash
#
gfortran -c -g sparsekit_prb11.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparsekit_prb11.f90"
  exit
fi
rm compiler.txt
#
gfortran sparsekit_prb11.o -L$HOME/lib/$ARCH -lsparsekit
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparsekit_prb11.o"
  exit
fi
rm sparsekit_prb11.o
#
mv a.out sparsekit_prb11
./sparsekit_prb11 > sparsekit_prb11_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparsekit_prb11"
  exit
fi
rm sparsekit_prb11
#
echo "Program output written to sparsekit_prb11_output.txt"
