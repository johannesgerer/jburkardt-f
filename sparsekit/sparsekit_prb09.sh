#!/bin/bash
#
gfortran -c -g sparsekit_prb09.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparsekit_prb09.f90"
  exit
fi
rm compiler.txt
#
gfortran sparsekit_prb09.o -L$HOME/lib/$ARCH -lsparsekit
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparsekit_prb09.o"
  exit
fi
rm sparsekit_prb09.o
#
mv a.out sparsekit_prb09
./sparsekit_prb09 > sparsekit_prb09_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparsekit_prb09"
  exit
fi
rm sparsekit_prb09
#
echo "Program output written to sparsekit_prb09_output.txt"
