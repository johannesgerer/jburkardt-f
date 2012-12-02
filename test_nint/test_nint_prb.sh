#!/bin/bash
#
gfortran -c -g test_nint_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_nint_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran test_nint_prb.o -L$HOME/lib/$ARCH -ltest_nint
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_nint_prb.o"
  exit
fi
rm test_nint_prb.o
#
mv a.out test_nint_prb
./test_nint_prb > test_nint_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_nint_prb"
  exit
fi
rm test_nint_prb
#
echo "Test program output written to test_nint_prb_output.txt."
