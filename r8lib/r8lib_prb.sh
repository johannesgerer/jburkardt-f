#!/bin/bash
#
gfortran -c -g r8lib_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling r8lib_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran r8lib_prb.o -L$HOME/lib/$ARCH -lr8lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading r8lib_prb.o"
  exit
fi
rm r8lib_prb.o
#
mv a.out r8lib_prb
./r8lib_prb > r8lib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running r8lib_prb"
  exit
fi
rm r8lib_prb
#
echo "Test program output written to r8lib_prb_output.txt."
