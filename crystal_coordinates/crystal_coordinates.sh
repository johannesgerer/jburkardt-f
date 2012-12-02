#!/bin/bash
#
gfortran -c -g crystal_coordinates.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling crystal_coordinates.f90"
  exit
fi
rm compiler.txt
#
gfortran crystal_coordinates.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading crystal_coordinates.o"
  exit
fi
rm crystal_coordinates.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/crystal_coordinates
#
echo "Program installed as ~/bin/$ARCH/crystal_coordinates"
