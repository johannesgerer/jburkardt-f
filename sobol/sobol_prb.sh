#!/bin/bash
#
gfortran -c -g sobol_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sobol_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sobol_prb.o -L$HOME/lib/$ARCH -lsobol
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sobol_prb.o"
  exit
fi
rm sobol_prb.o
#
mv a.out sobol_prb
./sobol_prb > sobol_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sobol_prb"
  exit
fi
rm sobol_prb
#
echo "Test program output written to sobol_prb_output.txt."
