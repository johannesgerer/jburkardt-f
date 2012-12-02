#!/bin/bash
#
g95 -c -g constant_type.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling constant_type.f90"
  exit
fi
rm compiler.txt
#
g95 constant_type.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading constant_type.o"
  exit
fi
rm constant_type.o
#
mv a.out constant_type_g95
./constant_type_g95 > constant_type_g95_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running constant_type_g95"
  exit
fi
rm constant_type_g95
#
echo "The constant_type_g95 test problem has been executed."
