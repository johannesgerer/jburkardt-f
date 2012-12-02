#!/bin/bash
#
gfortran -c -g toms661_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms661_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran toms661_prb.o -L$HOME/lib/$ARCH -ltoms661
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms661_prb.o"
  exit
fi
rm toms661_prb.o
#
mv a.out toms661_prb
./toms661_prb > toms661_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms661_prb"
  exit
fi
rm toms661_prb
#
echo "Test program output written to toms661_prb_output.txt."
