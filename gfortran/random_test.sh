#!/bin/bash
#
gfortran -c -g random_test.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling random_test.f90"
  exit
fi
rm compiler.txt
#
gfortran random_test.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading random_test.o"
  exit
fi
rm random_test.o
#
mv a.out random_test
./random_test > random_test_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running random_test"
  exit
fi
rm random_test
#
echo "Program output written to random_test_output.txt"
