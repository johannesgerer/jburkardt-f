#!/bin/bash
#
gfortran -c -g jacobi_ge_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling jacobi_ge_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran jacobi_ge_prb.o -L$HOME/lib/$ARCH -ltemplates
if [ $? -ne 0 ]; then
  echo "Errors linking and loading jacobi_ge_prb.o"
  exit
fi
rm jacobi_ge_prb.o
#
mv a.out jacobi_ge_prb
./jacobi_ge_prb > jacobi_ge_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running jacobi_ge_prb"
  exit
fi
rm jacobi_ge_prb
#
echo "Test program output written to jacobi_ge_prb_output.txt."
