#!/bin/bash
#
gfortran -c -g maxmin_test.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling maxmin_test.f90"
  exit
fi
rm compiler.txt
#
gfortran maxmin_test.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading maxmin_test.o"
  exit
fi
rm maxmin_test.o
#
mv a.out maxmin_test
./maxmin_test > maxmin_test_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running maxmin_test"
  exit
fi
rm maxmin_test
#
echo "Program output written maxmin_test_output.txt"
