#!/bin/bash
#
gfortran -c -g subanagram.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors while compiling subanagram.f90"
  exit
fi
rm compiler.txt
#
gfortran subanagram.o
if [ $? -ne 0 ]; then
  echo "Errors while loading subanagram.o"
  exit
fi
rm subanagram.o
#
mv a.out ~/bin/$ARCH/subanagram
#
echo "Executable installed as ~/bin/$ARCH/subanagram"
