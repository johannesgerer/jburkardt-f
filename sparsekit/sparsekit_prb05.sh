#!/bin/bash
#
gfortran -c -g sparsekit_prb05.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparsekit_prb05.f90"
  exit
fi
rm compiler.txt
#
gfortran sparsekit_prb05.o -L$HOME/lib/$ARCH -lsparsekit
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparsekit_prb05.o"
  exit
fi
rm sparsekit_prb05.o
#
mv a.out sparsekit_prb05
./sparsekit_prb05 > sparsekit_prb05_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparsekit_prb05"
  exit
fi
rm sparsekit_prb05
#
echo "Program output written to sparsekit_prb05_output.txt"
