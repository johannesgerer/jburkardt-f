#!/bin/bash
#
gfortran -c -g sparsekit_prb13.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparsekit_prb13.f90"
  exit
fi
rm compiler.txt
#
gfortran sparsekit_prb13.o -L$HOME/lib/$ARCH -lsparsekit
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparsekit_prb13.o"
  exit
fi
rm sparsekit_prb13.o
#
mv a.out sparsekit_prb13
./sparsekit_prb13 > sparsekit_prb13_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparsekit_prb13"
  exit
fi
rm sparsekit_prb13
#
echo "Program output written to sparsekit_prb13_output.txt"
