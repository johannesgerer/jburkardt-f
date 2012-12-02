#!/bin/bash
#
gfortran -c -g test_min_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_min_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran test_min_prb.o -L$HOME/lib/$ARCH -ltest_min
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_min_prb.o"
  exit
fi
rm test_min_prb.o
#
mv a.out test_min_prb
./test_min_prb > test_min_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_min_prb"
  exit
fi
rm test_min_prb
#
echo "Test program output written to test_min_prb_output.txt."
