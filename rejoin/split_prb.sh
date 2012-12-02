#!/bin/bash
#
gfortran -c -g split_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling split_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran split_prb.o -L$HOME/lib/$ARCH -lrejoin
if [ $? -ne 0 ]; then
  echo "Errors linking and loading split_prb.o"
  exit
fi
rm split_prb.o
#
mv a.out split_prb
./split_prb > split_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running split_prb"
  exit
fi
rm split_prb
#
echo "Test program output written to split_prb_output.txt."
