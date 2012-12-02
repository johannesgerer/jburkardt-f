#!/bin/bash
#
gfortran -c -g triangulation_rcm.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangulation_rcm.f90"
  exit
fi
rm compiler.txt
#
gfortran triangulation_rcm.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangulation_rcm.o"
  exit
fi
rm triangulation_rcm.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/triangulation_rcm
#
echo "Executable installed as ~/bin/$ARCH/triangulation_rcm"
