#!/bin/bash
#
gfortran -c -g nintlib_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling nintlib_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran nintlib_prb.o -L$HOME/lib/$ARCH -lnintlib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading nintlib_prb.o"
  exit
fi
rm nintlib_prb.o
#
mv a.out nintlib_prb
./nintlib_prb > nintlib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running nintlib_prb"
  exit
fi
rm nintlib_prb
#
echo "Test program output written to nintlib_prb_output.txt."
