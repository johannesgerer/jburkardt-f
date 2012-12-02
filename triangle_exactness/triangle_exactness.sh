#!/bin/bash
#
gfortran -c -g triangle_exactness.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_exactness.f90"
  exit
fi
rm compiler.txt
#
gfortran triangle_exactness.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_exactness.o"
  exit
fi
rm triangle_exactness.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/triangle_exactness
#
echo "Executable installed as ~/bin/$ARCH/triangle_exactness"
