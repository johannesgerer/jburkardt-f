#!/bin/bash
#
gfortran -c -g test_ls_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_ls_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran test_ls_prb.o -L$HOME/lib/$ARCH -ltest_ls -lr8lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_ls_prb.o"
  exit
fi
rm test_ls_prb.o
#
mv a.out test_ls_prb
./test_ls_prb > test_ls_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_ls_prb"
  exit
fi
rm test_ls_prb
#
echo "Test program output written to test_ls_prb_output.txt."
