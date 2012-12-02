#!/bin/bash
#
gfortran -c -g matrix_exponential_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling matrix_exponential_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran matrix_exponential_prb.o -L$HOME/lib/$ARCH -lmatrix_exponential \
  -ltest_matrix_exponential -lr8lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading matrix_exponential_prb.o"
  exit
fi
rm matrix_exponential_prb.o
#
mv a.out matrix_exponential_prb
./matrix_exponential_prb > matrix_exponential_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running matrix_exponential_prb"
  exit
fi
rm matrix_exponential_prb
#
echo "Test program output written to matrix_exponential_prb_output.txt."
