#!/bin/bash
#
gfortran -c division.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling division.f90"
  exit
fi
rm compiler.txt
#
gfortran division.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading division.o"
  exit
fi
rm division.o
#
mv a.out division_g95
./division_g95 > division_g95_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running division_g95"
  exit
fi
rm division_g95
#
echo "The division_g95 test problem has been executed."
