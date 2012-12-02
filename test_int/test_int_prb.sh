#!/bin/bash
#
gfortran -c -g test_int_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_int_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran test_int_prb.o -L$HOME/lib/$ARCH -ltest_int
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_int_prb.o"
  exit
fi
rm test_int_prb.o
#
mv a.out test_int_prb
./test_int_prb > test_int_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_int_prb"
  exit
fi
rm test_int_prb
#
echo "Test program output written to test_int_prb_output.txt."
