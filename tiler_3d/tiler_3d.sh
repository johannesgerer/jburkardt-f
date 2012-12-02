#!/bin/bash
#
gfortran -c -g tiler_3d.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling tiler_3d.f90"
  exit
fi
rm compiler.txt
#
gfortran tiler_3d.o
if [ $? -ne 0 ]; then
  echo "Errors loading tiler_3d.o"
  exit
fi
#
rm tiler_3d.o
mv a.out ~/bin/$ARCH/tiler_3d
#
echo "Executable installed as ~/bin/$ARCH/tiler_3d"
