#!/bin/bash
#
gfortran -c -g test_ode_prb.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling test_ode_prb.f90"
  exit
fi
rm compiler.txt
#
gfortran test_ode_prb.o -L$HOME/lib/$ARCH -ltest_ode
if [ $? -ne 0 ]; then
  echo "Errors linking and loading test_ode_prb.o"
  exit
fi
rm test_ode_prb.o
#
mv a.out test_ode_prb
./test_ode_prb > test_ode_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running test_ode_prb"
  exit
fi
rm test_ode_prb
#
echo "Test program output written to test_ode_prb_output.txt."
