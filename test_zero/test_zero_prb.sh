#!/bin/bash
#
gfortran -c -g test_zero_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_zero_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran test_zero_prb.o -L$HOME/lib/$ARCH -ltest_zero
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_zero_prb.o"
  exit
fi
rm test_zero_prb.o
#
mv a.out test_zero_prb
./test_zero_prb > test_zero_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_zero_prb"
  exit
fi
rm test_zero_prb
#
echo "Test program output written to test_zero_prb_output.txt."
