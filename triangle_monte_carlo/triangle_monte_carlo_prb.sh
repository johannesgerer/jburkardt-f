#!/bin/bash
#
gfortran -c -g triangle_monte_carlo_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_monte_carlo_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran triangle_monte_carlo_prb.o -L$HOME/lib/$ARCH -ltriangle_monte_carlo
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_monte_carlo_prb.o"
  exit
fi
rm triangle_monte_carlo_prb.o
#
mv a.out triangle_monte_carlo_prb
./triangle_monte_carlo_prb > triangle_monte_carlo_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running triangle_monte_carlo_prb"
  exit
fi
rm triangle_monte_carlo_prb
#
echo "Test program output written to triangle_monte_carlo_prb_output.txt."
