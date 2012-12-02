#!/bin/bash
#
gfortran -c -g triangulation_quad.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_quad.f90"
  exit
fi
rm compiler.txt
#
gfortran triangulation_quad.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_quad.o"
  exit
fi
rm triangulation_quad.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/triangulation_quad
#
echo "Executable installed as ~/bin/$ARCH/triangulation_quad"
