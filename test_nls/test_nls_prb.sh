#!/bin/bash
#
gfortran -c -g test_nls_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_nls_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran test_nls_prb.o -L$HOME/lib/$ARCH -ltest_nls
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_nls_prb.o"
  exit
fi
rm test_nls_prb.o
#
mv a.out test_nls_prb
./test_nls_prb > test_nls_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_nls_prb"
  exit
fi
rm test_nls_prb
#
echo "Test program output written to test_nls_prb_output.txt."
