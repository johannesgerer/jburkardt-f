#!/bin/bash
#
gfortran -c -g recursive_sub_test.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling recursive_sub_test.f90"
  exit
fi
rm compiler.txt
#
gfortran recursive_sub_test.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading recursive_sub_test.o"
  exit
fi
rm recursive_sub_test.o
#
mv a.out recursive_sub_test
./recursive_sub_test > recursive_sub_test_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running recursive_sub_test"
  exit
fi
rm recursive_sub_test
#
echo "The recursive_sub_test test problem has been executed."
