#!/bin/bash
#
gfortran -c -g subset_sum_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling subset_sum_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran subset_sum_prb.o -L$HOME/lib/$ARCH -lsubset_sum
if [ $? -ne 0 ]; then
  echo "Errors linking and loading subset_sum_prb.o"
  exit
fi
rm subset_sum_prb.o
#
mv a.out subset_sum_prb
./subset_sum_prb > subset_sum_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running subset_sum_prb"
  exit
fi
rm subset_sum_prb
#
echo "Test program output written to subset_sum_prb_output.txt."
