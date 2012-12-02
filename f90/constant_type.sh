#!/bin/bash
#
gfortran -c -g constant_type.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling constant_type.f90"
  exit
fi
rm compiler.txt
#
gfortran constant_type.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading constant_type.o"
  exit
fi
rm constant_type.o
#
mv a.out constant_type
./constant_type > constant_type_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running constant_type"
  exit
fi
rm constant_type
#
echo "The constant_type test problem has been executed."
