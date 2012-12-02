#!/bin/bash
#
gfortran -c -g test_nonlin_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_nonlin_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran test_nonlin_prb.o -L$HOME/lib/$ARCH -ltest_nonlin
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_nonlin_prb.o"
  exit
fi
rm test_nonlin_prb.o
#
mv a.out test_nonlin_prb
./test_nonlin_prb > test_nonlin_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_nonlin_prb"
  exit
fi
rm test_nonlin_prb
#
echo "Test program output written to test_nonlin_prb_output.txt."
