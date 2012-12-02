#!/bin/bash
#
gfortran -c -g triangulation_orient.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_orient.f90"
  exit
fi
rm compiler.txt
#
gfortran triangulation_orient.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_orient.o"
  exit
fi
#
rm triangulation_orient.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/triangulation_orient
#
echo "Executable installed as ~/bin/$ARCH/triangulation_orient"
