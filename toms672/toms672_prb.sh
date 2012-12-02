#!/bin/bash
#
gfortran -c -g toms672_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling toms672_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran toms672_prb.o -L$HOME/lib/$ARCH -ltoms672
if [ $? -ne 0 ]; then
  echo "Errors linking and loading toms672_prb.o"
  exit
fi
rm toms672_prb.o
#
mv a.out toms672_prb
./toms672_prb < toms672_prb_input.txt > toms672_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms672_prb"
  exit
fi
rm toms672_prb
#
echo "Test program output written to toms672_prb_output.txt."
