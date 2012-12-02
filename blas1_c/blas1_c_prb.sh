#!/bin/bash
#
gfortran -c -g blas1_c_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling blas1_c_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran blas1_c_prb.o -L$HOME/lib/$ARCH -lblas1_c
if [ $? -ne 0 ]; then
  echo "Errors linking and loading blas1_c_prb.o"
  exit
fi
rm blas1_c_prb.o
#
mv a.out blas1_c_prb
./blas1_c_prb > blas1_c_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running blas1_c_prb"
  exit
fi
rm blas1_c_prb
#
echo "Test program output written to blas1_c_prb_output.txt."
