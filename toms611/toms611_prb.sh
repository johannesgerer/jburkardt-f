#!/bin/bash
#
gfortran -c -g toms611_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms611_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran toms611_prb.o -L$HOME/lib/$ARCH -ltoms611
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms611_prb.o"
  exit
fi
rm toms611_prb.o
#
mv a.out toms611_prb
./toms611_prb > toms611_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms611_prb"
  exit
fi
rm toms611_prb
#
echo "Test program output written to toms611_prb_output.txt."
