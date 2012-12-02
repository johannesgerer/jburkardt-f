#!/bin/bash
#
gfortran -c -g cost_test.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cost_test.f90"
  exit
fi
rm compiler.txt
#
gfortran cost_test.o -L$HOME/lib/$ARCH -lfftpack5
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cost_test.o"
  exit
fi
rm cost_test.o
#
mv a.out cost_test
./cost_test > cost_test_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running cost_test"
  exit
fi
rm cost_test
#
echo "Test program output written to cost_test_output.txt."
