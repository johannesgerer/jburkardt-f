#!/bin/bash
#
gfortran -c -g digits.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling digits.f90"
  exit
fi
rm compiler.txt
#
gfortran digits.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading digits.o"
  exit
fi
rm digits.o
#
mv a.out digits
./digits > digits_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running digits"
  exit
fi
rm digits
#
echo "Program output was written to digits_output.txt"
