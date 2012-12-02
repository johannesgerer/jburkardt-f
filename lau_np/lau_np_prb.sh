#!/bin/bash
#
gfortran -c -g lau_np_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lau_np_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran lau_np_prb.o -L$HOME/lib/$ARCH -llau_np
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lau_np_prb.o"
  exit
fi
rm lau_np_prb.o
#
mv a.out lau_np_prb
./lau_np_prb > lau_np_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lau_np_prb"
  exit
fi
rm lau_np_prb
#
echo "Test program output written to lau_np_prb_output.txt."
