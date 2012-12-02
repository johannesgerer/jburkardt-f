#!/bin/bash
#
gfortran -c -g cvt_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran cvt_prb.o -L$HOME/lib/$ARCH -lcvt
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cvt_prb.o"
  exit
fi
rm cvt_prb.o
#
mv a.out cvt_prb
./cvt_prb > cvt_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cvt_prb"
  exit
fi
rm cvt_prb
#
echo "Test program output written to cvt_prb_output.txt."
