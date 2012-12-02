#!/bin/bash
#
gfortran -c -g jacobi_polynomial_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling jacobi_polynomial_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran jacobi_polynomial_prb.o -L$HOME/lib/$ARCH -ljacobi_polynomial
if [ $? -ne 0 ]; then
  echo "Errors linking and loading jacobi_polynomial_prb.o"
  exit
fi
rm jacobi_polynomial_prb.o
#
mv a.out jacobi_polynomial_prb
./jacobi_polynomial_prb > jacobi_polynomial_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running jacobi_polynomial_prb"
  exit
fi
rm jacobi_polynomial_prb
#
echo "Program output written to jacobi_polynomial_prb_output.txt"
