#!/bin/bash
#
gfortran -c -g triangle_properties.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling triangle_properties.f90"
  exit
fi
rm compiler.txt
#
gfortran triangle_properties.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading triangle_properties.o"
  exit
fi
rm triangle_properties.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/triangle_properties
#
echo "Executable installed as ~/bin/$ARCH/triangle_properties"
