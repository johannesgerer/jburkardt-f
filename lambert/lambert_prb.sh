#!/bin/bash
#
gfortran -c -g lambert_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lambert_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran lambert_prb.o -L$HOME/lib/$ARCH -llambert
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lambert_prb.o"
  exit
fi
rm lambert_prb.o
#
mv a.out lambert_prb
./lambert_prb > lambert_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running lambert_prb"
  exit
fi
rm lambert_prb
#
echo "Test program output written to lambert_prb_output.txt."
