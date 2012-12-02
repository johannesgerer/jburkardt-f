#!/bin/bash
#
gfortran -c -g cvt_fixed_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt_fixed_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran cvt_fixed_prb.o -L$HOME/lib/$ARCH -lcvt_fixed -lpbmlib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cvt_fixed_prb.o"
  exit
fi
rm cvt_fixed_prb.o
#
mv a.out cvt_fixed_prb
./cvt_fixed_prb > cvt_fixed_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cvt_fixed_prb"
  exit
fi
rm cvt_fixed_prb
#
echo "Test program output written to cvt_fixed_prb_output.txt."
