#!/bin/bash
#
gfortran -c -g cvt_weight_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt_weight_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran cvt_weight_prb.o -L$HOME/lib/$ARCH -lcvt_weight
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cvt_weight_prb.o"
  exit
fi
rm cvt_weight_prb.o
#
mv a.out cvt_weight_prb
./cvt_weight_prb > cvt_weight_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cvt_weight_prb"
  exit
fi
rm cvt_weight_prb
#
echo "Test program output written to cvt_weight_prb_output.txt."
