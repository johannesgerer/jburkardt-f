#!/bin/bash
#
gfortran -c -g -I/usr/local/dislin/gfortran dislin_ex07.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dislin_ex07.f90"
  exit
fi
rm compiler.txt
#
gfortran dislin_ex07.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dislin_ex07.o"
  exit
fi
rm dislin_ex07.o
#
mv a.out dislin_ex07
./dislin_ex07 > dislin_ex07_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dislin_ex07"
  exit
fi
rm dislin_ex07
#
echo "Program output written to dislin_ex07_output.txt"
