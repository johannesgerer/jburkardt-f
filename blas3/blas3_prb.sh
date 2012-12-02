#!/bin/bash
#
gfortran -c -g blas3_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas3_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran blas3_prb.o -L$HOME/lib/$ARCH -lblas3 -lblas2 -lblas1
if [ $? -ne 0 ]; then
  echo "Errors linking and loading blas3_prb.o"
  exit
fi
rm blas3_prb.o
#
mv a.out blas3_prb
./blas3_prb > blas3_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running blas3_prb"
  exit
fi
rm blas3_prb
#
echo "Test program output written to blas3_prb_output.txt."
