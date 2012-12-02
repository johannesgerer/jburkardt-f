#!/bin/bash
#
gfortran -c -g sparsekit_prb03.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparsekit_prb03.f90"
  exit
fi
rm compiler.txt
#
gfortran sparsekit_prb03.o -L$HOME/lib/$ARCH -lsparsekit
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparsekit_prb03.o"
  exit
fi
rm sparsekit_prb03.o
#
mv a.out sparsekit_prb03
./sparsekit_prb03 > sparsekit_prb03_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sparsekit_prb03"
  exit
fi
rm sparsekit_prb03
#
echo "Program output written to sparsekit_prb03_output.txt"
