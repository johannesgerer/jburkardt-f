#!/bin/bash
#
gfortran -c -g toms790_prb2.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms790_prb2.f90"
  exit
fi
rm compiler.txt
#
gfortran toms790_prb2.o -L$HOME/lib/$ARCH -ltoms790
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms790_prb2.o"
  exit
fi
rm toms790_prb2.o
#
mv a.out toms790_prb2
./toms790_prb2 > toms790_prb2_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms790_prb2"
  exit
fi
rm toms790_prb2
#
echo "Test results written to toms790_prb2_output.txt."
