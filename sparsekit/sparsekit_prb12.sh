#!/bin/bash
#
gfortran -c -g sparsekit_prb12.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparsekit_prb12.f90"
  exit
fi
rm compiler.txt
#
gfortran sparsekit_prb12.o -L$HOME/lib/$ARCH -lsparsekit
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparsekit_prb12.o"
  exit
fi
rm sparsekit_prb12.o
#
mv a.out sparsekit_prb12
./sparsekit_prb12 > sparsekit_prb12_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparsekit_prb12"
  exit
fi
rm sparsekit_prb12
#
echo "Program output written to sparsekit_prb12_output.txt"
