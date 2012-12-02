#!/bin/bash
#
gfortran -c -g genin_two.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Error while compiling genin_two.f90"
  exit
fi
rm compiler.txt
#
gfortran genin_two.o
if [ $? -ne 0 ]; then
  echo "Error while loading genin_two.o"
  exit
fi
rm genin_two.o
#
mv a.out ~/bin/$ARCH/genin_two
#
echo "Executable installed as ~/bin/$ARCH/genin_two"
