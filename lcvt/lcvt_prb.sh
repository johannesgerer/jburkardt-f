#!/bin/bash
#
gfortran -c -g lcvt_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lcvt_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran lcvt_prb.o -L$HOME/lib/$ARCH -llcvt
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lcvt_prb.o"
  exit
fi
rm lcvt_prb.o
#
mv a.out lcvt_prb
./lcvt_prb > lcvt_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lcvt_prb"
  exit
fi
rm lcvt_prb
#
echo "Test program output written to lcvt_prb_output.txt."
