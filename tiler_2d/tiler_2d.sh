#!/bin/bash
#
gfortran -c -g tiler_2d.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tiler_2d.f90"
  exit
fi
rm compiler.txt
#
gfortran tiler_2d.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading tiler_2d.o"
  exit
fi
#
rm tiler_2d.o
mv a.out ~/bin/$ARCH/tiler_2d
#
echo "Executable installed as ~/bin/$ARCH/tiler_2d"
