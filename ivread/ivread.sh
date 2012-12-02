#!/bin/bash
#
gfortran -c -g ivread.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling ivread.f90"
  exit
fi
rm compiler.txt
#
gfortran ivread.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading ivread.o"
  exit
fi
#
rm ivread.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/ivread
#
echo "Executable installed as ~/bin/$ARCH/ivread"
