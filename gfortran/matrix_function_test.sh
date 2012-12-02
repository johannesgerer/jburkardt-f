#!/bin/bash
#
gfortran -c -g matrix_function_test.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling matrix_function_test.f90"
  exit
fi
rm compiler.txt
#
gfortran matrix_function_test.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading matrix_function_test.o"
  exit
fi
rm matrix_function_test.o
#
mv a.out matrix_function_test
./matrix_function_test > matrix_function_test_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running matrix_function_test"
  exit
fi
rm matrix_function_test
#
echo "Program output written to matrix_function_test_output.txt"
