#!/bin/bash
#
gfortran -c -g sparsekit_prb04.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparsekit_prb04.f90"
  exit
fi
rm compiler.txt
#
gfortran sparsekit_prb04.o -L$HOME/lib/$ARCH -lsparsekit
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparsekit_prb04.o"
  exit
fi
rm sparsekit_prb04.o
#
mv a.out sparsekit_prb04
./sparsekit_prb04 > sparsekit_prb04_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparsekit_prb04"
  exit
fi
rm sparsekit_prb04
#
echo "Program output written to sparsekit_prb04_output.txt"
