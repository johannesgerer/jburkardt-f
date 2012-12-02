#!/bin/bash
#
gfortran -c -g intlib_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling intlib_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran intlib_prb.o -L$HOME/lib/$ARCH -lintlib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading intlib_prb.o"
  exit
fi
rm intlib_prb.o
#
mv a.out intlib_prb
./intlib_prb > intlib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running intlib_prb"
  exit
fi
rm intlib_prb
#
echo "Test program output written to intlib_prb_output.txt."
