#!/bin/bash
#
gfortran -c -g -I/usr/local/dislin/gfortran dislin_ex05.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dislin_ex05.f90"
  exit
fi
rm compiler.txt
#
gfortran dislin_ex05.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dislin_ex05.o"
  exit
fi
rm dislin_ex05.o
#
mv a.out dislin_ex05
./dislin_ex05 > dislin_ex05_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dislin_ex05"
  exit
fi
rm dislin_ex05
#
echo "Program output written to dislin_ex05_output.txt"
