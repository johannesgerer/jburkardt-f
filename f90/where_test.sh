#!/bin/bash
#
gfortran -c -g where_test.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling where_test.f90"
  exit
fi
rm compiler.txt
#
gfortran where_test.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading where_test.o"
  exit
fi
rm where_test.o
#
mv a.out where_test
./where_test > where_test_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running where_test"
  exit
fi
rm where_test
#
echo "The where_test test problem has been executed."
