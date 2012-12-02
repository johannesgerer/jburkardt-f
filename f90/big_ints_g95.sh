#!/bin/bash
#
gfortran -c -g big_ints.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling big_ints.f90"
  exit
fi
rm compiler.txt
#
gfortran big_ints.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading big_ints.o"
  exit
fi
rm big_ints.o
#
mv a.out big_ints_g95
./big_ints_g95 > big_ints_g95_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running big_ints_g95"
  exit
fi
rm big_ints_g95
#
echo "The big_ints_g95 test problem has been executed."
