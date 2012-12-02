#!/bin/bash
#
gfortran -c circle_test.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Error compiling circle_test.f90"
  exit
fi
rm compiler.txt
#
gfortran circle_test.o
if [ $? -ne 0 ]; then
  echo "Error loading circle_test.o"
  exit
fi
rm circle_test.o
#
mv a.out ~/bin/$ARCH/circle_test
#
echo "Executable installed as ~/bin/$ARCH/circle_test"
