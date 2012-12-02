#!/bin/bash
#
gfortran -c -g crystal.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling crystal.f90"
  exit
fi
rm compiler.txt
#
gfortran crystal.o -L$HOME/lib/$ARCH -ltoms611
if [ $? -ne 0 ]; then
  echo "Errors while linking and loading crystal.o"
  exit
fi
rm crystal.o
#
mv a.out ~/bin/$ARCH/crystal
#
echo "Executable installed as ~/bin/$ARCH/crystal"
