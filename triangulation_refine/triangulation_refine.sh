#!/bin/bash
#
gfortran -c -g triangulation_refine.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_refine.f90"
  exit
fi
rm compiler.txt
#
gfortran triangulation_refine.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_refine.o"
  exit
fi
#
rm triangulation_refine.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/triangulation_refine
#
echo "Executable installed as ~/bin/$ARCH/triangulation_refine"
