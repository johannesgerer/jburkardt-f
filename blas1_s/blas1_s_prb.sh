#!/bin/bash
#
gfortran -c -g blas1_s_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas1_s_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran blas1_s_prb.o -L$HOME/lib/$ARCH -lblas1_s
if [ $? -ne 0 ]; then
  echo "Errors linking and loading blas1_s_prb.o"
  exit
fi
rm blas1_s_prb.o
#
mv a.out blas1_s_prb
./blas1_s_prb > blas1_s_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running blas1_s_prb"
  exit
fi
rm blas1_s_prb
#
echo "Test program output written to blas1_s_prb_output.txt."
