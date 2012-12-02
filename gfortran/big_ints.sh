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
mv a.out big_ints
./big_ints > big_ints_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running big_ints"
  exit
fi
rm big_ints
#
echo "Program output written to big_ints_output.txt"
