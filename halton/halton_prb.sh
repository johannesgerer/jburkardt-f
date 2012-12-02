#!/bin/bash
#
gfortran -c -g halton_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling halton_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran halton_prb.o -L$HOME/lib/$ARCH -lhalton
if [ $? -ne 0 ]; then
  echo "Errors linking and loading halton_prb.o"
  exit
fi
rm halton_prb.o
#
mv a.out halton_prb
./halton_prb > halton_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running halton_prb"
  exit
fi
rm halton_prb
#
echo "Test program output written to halton_prb_output.txt."
