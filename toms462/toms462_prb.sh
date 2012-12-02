#!/bin/bash
#
gfortran -c -g toms462_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms462_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran toms462_prb.o -L$HOME/lib/$ARCH -ltoms462
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms462_prb.o"
  exit
fi
rm toms462_prb.o
#
mv a.out toms462_prb
./toms462_prb > toms462_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms462_prb"
  exit
fi
rm toms462_prb
#
echo "Test program output written to toms462_prb_output.txt."
