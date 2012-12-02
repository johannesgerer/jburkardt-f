#!/bin/bash
#
gfortran -c -g regression_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling regression_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran regression_prb.o -L$HOME/lib/$ARCH -lregression
if [ $? -ne 0 ]; then
  echo "Errors linking and loading regression_prb.o"
  exit
fi
rm regression_prb.o
#
mv a.out regression_prb
./regression_prb > regression_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running regression_prb"
  exit
fi
rm regression_prb
#
echo "Test program output written to regression_prb_output.txt."
