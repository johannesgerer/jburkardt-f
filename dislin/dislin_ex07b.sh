#!/bin/bash
#
gfortran -c -g -I/usr/local/dislin/gfortran dislin_ex07b.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling dislin_ex07b.f90"
  exit
fi
rm compiler.txt
#
gfortran dislin_ex07b.o -L/usr/local/dislin -ldislin -L/opt/local/lib -lXm
if [ $? -ne 0 ]; then
  echo "Errors linking and loading dislin_ex07b.o"
  exit
fi
rm dislin_ex07b.o
#
mv a.out dislin_ex07b
./dislin_ex07b > dislin_ex07b_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running dislin_ex07b"
  exit
fi
rm dislin_ex07b
#
echo "Program output written to dislin_ex07b_output.txt"
