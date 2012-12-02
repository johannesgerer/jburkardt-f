#!/bin/bash
#
gfortran -c -g piecewise_linear_product_integral_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling piecewise_linear_product_integral_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran piecewise_linear_product_integral_prb.o -L$HOME/lib/$ARCH -lpiecewise_linear_product_integral
if [ $? -ne 0 ]; then
  echo "Errors linking and loading piecewise_linear_product_integral_prb.o"
  exit
fi
rm piecewise_linear_product_integral_prb.o
#
mv a.out piecewise_linear_product_integral_prb
./piecewise_linear_product_integral_prb > piecewise_linear_product_integral_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running piecewise_linear_product_integral_prb"
  exit
fi
rm piecewise_linear_product_integral_prb
#
echo "Test program output written to piecewise_linear_product_integral_prb_output.txt."
