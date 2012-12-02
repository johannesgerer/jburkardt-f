#!/bin/bash
#
gfortran -c -g zero_one.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling zero_one.f90"
  exit
fi
rm compiler.txt
#
gfortran zero_one.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading zero_one.o"
  exit
fi
rm zero_one.o
#
mv a.out zero_one
#
./zero_one > zero_one_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running zero_one"
  exit
fi
rm zero_one
#
echo "Program output written to zero_one_output.txt"
