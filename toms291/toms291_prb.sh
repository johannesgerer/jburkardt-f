#!/bin/bash
#
gfortran -c -g toms291_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms291_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran toms291_prb.o -L$HOME/lib/$ARCH -ltoms291
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms291_prb.o"
  exit
fi
rm toms291_prb.o
#
mv a.out toms291_prb
./toms291_prb > toms291_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms291_prb"
  exit
fi
rm toms291_prb
#
echo "Test program output written to toms291_prb_output.txt."
