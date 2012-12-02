#!/bin/bash
#
gfortran -c -g game31.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling game31.f90"
  exit
fi
rm compiler.txt
#
gfortran game31.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading game31.o"
  exit
fi
rm game31.o
#
mv a.out ~/bin/$ARCH/game31
#
echo "Executable installed as ~/bin/$ARCH/game31"
