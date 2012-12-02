#!/bin/bash
#
gfortran -c -g exponent_format_overflow.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling exponent_format_overflow.f90"
  exit
fi
rm compiler.txt
#
gfortran exponent_format_overflow.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading exponent_format_overflow.o"
  exit
fi
rm exponent_format_overflow.o
#
mv a.out exponent_format_overflow
./exponent_format_overflow > exponent_format_overflow_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running exponent_format_overflow"
  exit
fi
rm exponent_format_overflow
#
echo "Program output written to exponent_format_overflow_output.txt"
