#!/bin/bash
#
gfortran -c -g toms726_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms726_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran toms726_prb.o -L$HOME/lib/$ARCH -ltoms726
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms726_prb.o"
  exit
fi
rm toms726_prb.o
#
mv a.out toms726_prb
./toms726_prb > toms726_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms726_prb"
  exit
fi
rm toms726_prb
#
echo "Test program output written to toms726_prb_output.txt."
