#!/bin/bash
#
gfortran -c -g interp_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling interp_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran interp_prb.o -L$HOME/lib/$ARCH -linterp
if [ $? -ne 0 ]; then
  echo "Errors linking and loading interp_prb.o"
  exit
fi
rm interp_prb.o
#
mv a.out interp_prb
./interp_prb > interp_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running interp_prb"
  exit
fi
rm interp_prb
#
echo "Test program output written to interp_prb_output.txt."
