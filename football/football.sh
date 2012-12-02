#!/bin/bash
#
gfortran -c football.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling football.f90"
  exit
fi
rm compiler.txt
#
gfortran football.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading football.o"
  exit
fi
rm football.o
#
mv a.out ~/bin/$ARCH/football
#
echo "Program installed as ~/bin/$ARCH/football"
