#!/bin/bash
#
gfortran -c -g sparsekit_prb06.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparsekit_prb06.f90"
  exit
fi
rm compiler.txt
#
gfortran sparsekit_prb06.o -L$HOME/lib/$ARCH -lsparsekit
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparsekit_prb06.o"
  exit
fi
rm sparsekit_prb06.o
#
mv a.out sparsekit_prb06
./sparsekit_prb06 > sparsekit_prb06_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparsekit_prb06"
  exit
fi
rm sparsekit_prb06
#
echo "Program output written to sparsekit_prb06_output.txt"
