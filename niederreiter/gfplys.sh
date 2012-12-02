#!/bin/bash
#
gfortran -c -g gfplys.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Error while compiling gfplys.f90"
  exit
fi
rm compiler.txt
#
gfortran gfplys.o
if [ $? -ne 0 ]; then
  echo "Error while loading gfplys.o"
  exit
fi
rm gfplys.o
#
mv a.out ~/bin/$ARCH/gfplys
#
echo "Executable installed as ~/bin/$ARCH/gfplys"
