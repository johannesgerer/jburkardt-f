#!/bin/bash
#
gfortran -c -g cvt_size_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt_size_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran cvt_size_prb.o -L$HOME/lib/$ARCH -lcvt_size
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cvt_size_prb.o"
  exit
fi
rm cvt_size_prb.o
#
mv a.out cvt_size_prb
./cvt_size_prb > cvt_size_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cvt_size_prb"
  exit
fi
rm cvt_size_prb
#
echo "Test program output written to cvt_size_prb_output.txt."
