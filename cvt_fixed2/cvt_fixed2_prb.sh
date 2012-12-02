#!/bin/bash
#
gfortran -c -g cvt_fixed2_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt_fixed2_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran cvt_fixed2_prb.o -L$HOME/lib/$ARCH -lcvt_fixed2
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cvt_fixed2_prb.o"
  exit
fi
rm cvt_fixed2_prb.o
#
mv a.out cvt_fixed2_prb
./cvt_fixed2_prb > cvt_fixed2_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cvt_fixed2_prb"
  exit
fi
rm cvt_fixed2_prb
#
echo "Test program output written to cvt_fixed2_prb_output.txt."
