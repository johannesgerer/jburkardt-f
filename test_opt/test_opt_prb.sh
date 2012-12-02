#!/bin/bash
#
gfortran -c -g test_opt_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_opt_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran test_opt_prb.o -L$HOME/lib/$ARCH -ltest_opt
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_opt_prb.o"
  exit
fi
rm test_opt_prb.o
#
mv a.out test_opt_prb
./test_opt_prb > test_opt_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_opt_prb"
  exit
fi
rm test_opt_prb
#
echo "Test program output written to test_opt_prb_output.txt."
