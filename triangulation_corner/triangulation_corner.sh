#!/bin/bash
#
gfortran -c -g triangulation_corner.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_corner.f90"
  exit
fi
rm compiler.txt
#
gfortran triangulation_corner.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_corner.o"
  exit
fi
#
rm triangulation_corner.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/triangulation_corner
#
echo "Executable installed as ~/bin/$ARCH/triangulation_corner"
