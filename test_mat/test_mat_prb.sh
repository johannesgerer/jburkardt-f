#!/bin/bash
#
gfortran -c -g test_mat_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_mat_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran test_mat_prb.o -L$HOME/lib/$ARCH -ltest_mat
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_mat_prb.o"
  exit
fi
rm test_mat_prb.o
#
mv a.out test_mat_prb
./test_mat_prb > test_mat_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_mat_prb"
  exit
fi
rm test_mat_prb
#
echo "Test program output written to test_mat_prb_output.txt."
