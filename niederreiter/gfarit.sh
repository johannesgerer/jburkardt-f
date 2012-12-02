#!/bin/bash
#
gfortran -c -g gfarit.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Error while compiling gfarit.f90"
  exit
fi
rm compiler.txt
#
gfortran gfarit.o
if [ $? -ne 0 ]; then
  echo "Error while loading gfarit.o"
  exit
fi
rm gfarit.o
#
mv a.out ~/bin/$ARCH/gfarit
#
echo "Executable installed as ~/bin/$ARCH/gfarit"
