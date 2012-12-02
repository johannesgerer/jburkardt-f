#!/bin/bash
#
gfortran -c -g toms453_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms453_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran toms453_prb.o -L$HOME/lib/$ARCH -ltoms453
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms453_prb.o"
  exit
fi
rm toms453_prb.o
#
mv a.out toms453_prb
./toms453_prb > toms453_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms453_prb"
  exit
fi
rm toms453_prb
#
echo "Test program output written to toms453_prb_output.txt."
