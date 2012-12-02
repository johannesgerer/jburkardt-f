#!/bin/bash
#
gfortran -c -g cvt_refine_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt_refine_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran cvt_refine_prb.o -L$HOME/lib/$ARCH -lcvt_refine
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cvt_refine_prb.o"
  exit
fi
rm cvt_refine_prb.o
#
mv a.out cvt_refine_prb
./cvt_refine_prb > cvt_refine_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cvt_refine_prb"
  exit
fi
rm cvt_refine_prb
#
echo "Test program output written to cvt_refine_prb_output.txt."
