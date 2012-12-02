#!/bin/bash
#
gfortran -c -g sparse_grid_cc_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sparse_grid_cc_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran sparse_grid_cc_dataset.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sparse_grid_cc_dataset.o"
  exit
fi
rm sparse_grid_cc_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/sparse_grid_cc_dataset
#
echo "Executable installed as ~/bin/$ARCH/sparse_grid_cc_dataset"
