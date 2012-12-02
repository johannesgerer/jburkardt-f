#!/bin/bash
#
gfortran -c -g laguerre_polynomial_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling laguerre_polynomial_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran laguerre_polynomial_prb.o -L$HOME/lib/$ARCH -llaguerre_polynomial
if [ $? -ne 0 ]; then
  echo "Errors linking and loading laguerre_polynomial_prb.o"
  exit
fi
rm laguerre_polynomial_prb.o
#
mv a.out laguerre_polynomial_prb
./laguerre_polynomial_prb > laguerre_polynomial_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running laguerre_polynomial_prb"
  exit
fi
rm laguerre_polynomial_prb
#
echo "Test program output written to laguerre_polynomial_prb_output.txt."
