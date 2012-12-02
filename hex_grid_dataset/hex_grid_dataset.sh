#!/bin/bash
#
gfortran -c -g hex_grid_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hex_grid_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran hex_grid_dataset.o -L$HOME/lib/$ARCH -lhex_grid
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hex_grid_dataset.o"
  exit
fi
rm hex_grid_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/hex_grid_dataset
#
echo "Executable installed as ~/bin/$ARCH/hex_grid_dataset"
