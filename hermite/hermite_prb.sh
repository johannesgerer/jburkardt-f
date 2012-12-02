#!/bin/bash
#
gfortran -c -g hermite_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran hermite_prb.o -L$HOME/lib/$ARCH -lhermite
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hermite_prb.o"
  exit
fi
rm hermite_prb.o
#
mv a.out hermite_prb
./hermite_prb > hermite_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hermite_prb"
  exit
fi
rm hermite_prb
#
echo "Program output written to hermite_prb_output.txt"
