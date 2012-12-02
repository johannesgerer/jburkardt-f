#!/bin/bash
#
gfortran -c -g test_interp_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_interp_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran test_interp_prb.o -L$HOME/lib/$ARCH -ltest_interp -lr8lib
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_interp_prb.o"
  exit
fi
rm test_interp_prb.o
#
mv a.out test_interp_prb
./test_interp_prb > test_interp_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_interp_prb"
  exit
fi
rm test_interp_prb
#
echo "Test program output written to test_interp_prb_output.txt."
