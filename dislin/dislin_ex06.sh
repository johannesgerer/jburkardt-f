#!/bin/bash
#
gfortran -c -g -I/usr/local/dislin/gfortran dislin_ex06.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dislin_ex06.f90"
  exit
fi
rm compiler.txt
#
gfortran dislin_ex06.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dislin_ex06.o"
  exit
fi
rm dislin_ex06.o
#
mv a.out dislin_ex06
./dislin_ex06 > dislin_ex06_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dislin_ex06"
  exit
fi
rm dislin_ex06
#
echo "Program output written to dislin_ex06_output.txt"
