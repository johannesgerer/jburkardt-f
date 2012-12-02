#!/bin/bash
#
gfortran -c -g blas2_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas2_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran blas2_prb.o -L$HOME/lib/$ARCH -lblas2 -lblas1_s
if [ $? -ne 0 ]; then
  echo "Errors linking and loading blas2_prb.o"
  exit
fi
rm blas2_prb.o
#
mv a.out blas2_prb
./blas2_prb > blas2_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running blas2_prb"
  exit
fi
rm blas2_prb
#
echo "Test program output written to blas2_prb_output.txt."
