#!/bin/bash
#
gfortran -c -g normal_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling normal_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran normal_prb.o -L$HOME/lib/$ARCH -lnormal
if [ $? -ne 0 ]; then
  echo "Errors linking and loading normal_prb.o"
  exit
fi
rm normal_prb.o
#
mv a.out normal_prb
./normal_prb > normal_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running normal_prb"
  exit
fi
rm normal_prb
#
echo "Test program output written to normal_prb_output.txt."
