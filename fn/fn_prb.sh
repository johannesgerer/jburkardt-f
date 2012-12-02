#!/bin/bash
#
gfortran -c -g fn_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling fn_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran fn_prb.o -L$HOME/lib/$ARCH -lfn -ltest_values
if [ $? -ne 0 ]; then
  echo "Errors linking and loading fn_prb.o"
  exit
fi
rm fn_prb.o
#
mv a.out fn_prb
./fn_prb > fn_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running fn_prb"
  exit
fi
rm fn_prb
#
echo "Test results written to fn_prb_output.txt."
