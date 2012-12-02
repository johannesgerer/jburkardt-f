#!/bin/bash
#
gfortran -c -g chebyshev_polynomial_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling chebyshev_polynomial_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran chebyshev_polynomial_prb.o -L$HOME/lib/$ARCH -lchebyshev_polynomial
if [ $? -ne 0 ]; then
  echo "Errors linking and loading chebyshev_polynomial_prb.o"
  exit
fi
rm chebyshev_polynomial_prb.o
#
mv a.out chebyshev_polynomial_prb
./chebyshev_polynomial_prb > chebyshev_polynomial_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running chebyshev_polynomial_prb"
  exit
fi
rm chebyshev_polynomial_prb
#
echo "Program output written to chebyshev_polynomial_prb_output.txt"
