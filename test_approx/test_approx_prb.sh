#!/bin/bash
#
gfortran -c -g test_approx_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_approx_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran test_approx_prb.o -L$HOME/lib/$ARCH -ltest_approx -lspline
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_approx_prb.o"
  exit
fi
rm test_approx_prb.o
#
mv a.out test_approx_prb
./test_approx_prb > test_approx_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_approx_prb"
  exit
fi
rm test_approx_prb
#
echo "Test program output written to test_approx_prb_output.txt."
