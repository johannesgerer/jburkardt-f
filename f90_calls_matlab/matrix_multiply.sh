#!/bin/bash
#
gfortran -c -g matrix_multiply.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling matrix_multiply.f90"
  exit
fi
rm compiler.txt
#
gfortran matrix_multiply.o
if [ $? -ne 0 ]; then
  echo "Errors while loading matrix_multiply.o"
  exit
fi
rm matrix_multiply.o
#
mv a.out ~/bin/$ARCH/matrix_multiply
#
echo "A new version of matrix_multiply has been created."
