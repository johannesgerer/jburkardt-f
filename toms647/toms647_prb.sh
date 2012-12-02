#!/bin/bash
#
gfortran -c -g toms647_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms647_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran toms647_prb.o -L$HOME/lib/$ARCH -ltoms647
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms647_prb.o"
  exit
fi
rm toms647_prb.o
#
mv a.out toms647_prb
./toms647_prb > toms647_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms647_prb"
  exit
fi
rm toms647_prb
#
echo "Test program output written to toms647_prb_output.txt."
