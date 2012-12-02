#!/bin/bash
#
gfortran -c -g division.f90 >& compiler.txt
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
mv a.out division
./division > division_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running division"
  exit
fi
rm division
#
echo "The division test problem has been executed."
