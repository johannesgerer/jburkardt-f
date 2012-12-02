#!/bin/bash
#
gfortran -c -g sine_transform_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sine_transform_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran sine_transform_prb.o -L$HOME/lib/$ARCH -lsine_transform
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sine_transform_prb.o"
  exit
fi
rm sine_transform_prb.o
#
mv a.out sine_transform_prb
./sine_transform_prb > sine_transform_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sine_transform_prb"
  exit
fi
rm sine_transform_prb
#
echo "Test program output written to sine_transform_prb_output.txt."
