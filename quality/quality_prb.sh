#!/bin/bash
#
gfortran -c -g quality_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling quality_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran quality_prb.o -L$HOME/lib/$ARCH -lquality
if [ $? -ne 0 ]; then
  echo "Errors linking and loading quality_prb.o"
  exit
fi
rm quality_prb.o
#
mv a.out quality_prb
./quality_prb > quality_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running quality_prb"
  exit
fi
rm quality_prb
#
echo "Test program output written to quality_prb_output.txt."
