#!/bin/bash
#
gfortran -c -g shallow_water_1d.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling shallow_water_1d.f90"
  exit
fi
rm compiler.txt
#
gfortran shallow_water_1d.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading shallow_water_1d.o"
  exit
fi
rm shallow_water_1d.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/shallow_water_1d
#
echo "Program installed as ~/bin/$ARCH/shallow_water_1d"
