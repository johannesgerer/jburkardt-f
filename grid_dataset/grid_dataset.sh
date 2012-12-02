#!/bin/bash
#
gfortran -c -g grid_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling grid_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran grid_dataset.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading grid_dataset.o"
  exit
fi
rm grid_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/grid_dataset
#
echo "Executable installed as ~/bin/$ARCH/grid_dataset"
