#!/bin/bash
#
gfortran -c -g test_int_laguerre_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_int_laguerre_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran test_int_laguerre_prb.o -L$HOME/lib/$ARCH -ltest_int_laguerre
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_int_laguerre_prb.o"
  exit
fi
rm test_int_laguerre_prb.o
#
mv a.out test_int_laguerre_prb
./test_int_laguerre_prb > test_int_laguerre_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_int_laguerre_prb"
  exit
fi
rm test_int_laguerre_prb
#
echo "Test program output written to test_int_laguerre_prb_output.txt."
