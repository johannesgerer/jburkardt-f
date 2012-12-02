#!/bin/bash
#
gfortran -c -g recursive_fun_test.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling recursive_fun_test.f90"
  exit
fi
rm compiler.txt
#
gfortran recursive_fun_test.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading recursive_fun_test.o"
  exit
fi
rm recursive_fun_test.o
#
mv a.out recursive_fun_test
./recursive_fun_test > recursive_fun_test_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running recursive_fun_test"
  exit
fi
rm recursive_fun_test
#
echo "The recursive_fun_test test problem has been executed."
