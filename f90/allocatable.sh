#!/bin/bash
#
gfortran -c -g -std=f2003 allocatable.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling allocatable.f90"
  exit
fi
rm compiler.txt
#
gfortran -std=2003 allocatable.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading allocatable.o"
  exit
fi
rm allocatable.o
#
mv a.out allocatable
./allocatable > allocatable_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running allocatable"
  exit
fi
rm allocatable
#
echo "Program output written to allocatable_output.txt"
