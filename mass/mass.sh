#!/bin/bash
#
gfortran -c mass.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mass.f90"
  exit
fi
rm compiler.txt
#
gfortran mass.o -L$HOME/lib/$ARCH -llapack
if [ $? -ne 0 ]; then
  echo "Errors linking and loading mass.o"
  exit
fi
rm mass.o
#
mv a.out ~/bin/$ARCH/mass
#
echo "Program installed as ~/bin/$ARCH/mass"
