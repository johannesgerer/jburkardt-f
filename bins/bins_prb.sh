#!/bin/bash
#
gfortran -c -g bins_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bins_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran bins_prb.o -L$HOME/lib/$ARCH -lbins
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bins_prb.o"
  exit
fi
rm bins_prb.o
#
mv a.out bins_prb
./bins_prb > bins_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running bins_prb"
  exit
fi
rm bins_prb
#
echo "Test program output written to bins_prb_output.txt."
