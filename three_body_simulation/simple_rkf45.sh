#!/bin/bash
#
gfortran -c -g simple_rkf45.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling simple_rkf45.f90"
  exit
fi
rm compiler.txt
#
gfortran simple_rkf45.o -L$HOME/lib/$ARCH -lrkf45
if [ $? -ne 0 ]; then
  echo "Errors linking and loading simple_rkf45.o"
  exit
fi
rm simple_rkf45.o
#
mv a.out simple_rkf45
./simple_rkf45 > simple_rkf45_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running simple_rkf45"
  exit
fi
rm simple_rkf45
#
echo "Test program output written to simple_rkf45_output.txt."
