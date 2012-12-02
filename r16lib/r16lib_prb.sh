#!/bin/bash
#
gfortran -c -g r16lib_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling r16lib_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran r16lib_prb.o -L$HOME/lib/$ARCH -lr16lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading r16lib_prb.o"
  exit
fi
rm r16lib_prb.o
#
mv a.out r16lib_prb
./r16lib_prb > r16lib_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running r16lib_prb"
  exit
fi
rm r16lib_prb
#
echo "Test program output written to r16lib_prb_output.txt."
