#!/bin/bash
#
gfortran -c -g black_scholes_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling black_scholes_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran black_scholes_prb.o -L$HOME/lib/$ARCH -lblack_scholes
if [ $? -ne 0 ]; then
  echo "Errors linking and loading black_scholes_prb.o"
  exit
fi
rm black_scholes_prb.o
#
mv a.out black_scholes_prb
./black_scholes_prb > black_scholes_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running black_scholes_prb"
  exit
fi
rm black_scholes_prb
#
echo "Program output written to black_scholes_prb_output.txt"
