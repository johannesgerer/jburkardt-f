#!/bin/bash
#
gfortran -c -g hermite_cubic_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hermite_cubic_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran hermite_cubic_prb.o -L$HOME/lib/$ARCH -lhermite_cubic
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hermite_cubic_prb.o"
  exit
fi
rm hermite_cubic_prb.o
#
mv a.out hermite_cubic_prb
./hermite_cubic_prb > hermite_cubic_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running hermite_cubic_prb"
  exit
fi
rm hermite_cubic_prb
#
echo "Test program output written to hermite_cubic_prb_output.txt."
