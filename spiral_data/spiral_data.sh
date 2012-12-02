#!/bin/bash
#
gfortran -c -g spiral_data.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling spiral_data.f90"
  exit
fi
rm compiler.txt
#
gfortran spiral_data.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading spiral_data.o"
  exit
fi
#
rm spiral_data.o
mv a.out ~/bin/$ARCH/spiral_data
#
echo "Executable installed as ~/bin/ARCH/spiral_data"
