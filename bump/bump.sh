#!/bin/bash
#
gfortran -c -g bump.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling bump.f90"
  exit
fi
rm compiler.txt
#
gfortran bump.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading bump.o"
  exit
fi
rm bump.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/bump
#
echo "Program installed as ~/bin/$ARCH/bump"
