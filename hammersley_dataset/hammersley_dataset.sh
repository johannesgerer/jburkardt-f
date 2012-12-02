#!/bin/bash
#
gfortran -c -g hammersley_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling hammersley_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran hammersley_dataset.o -L$HOME/lib/$ARCH -lhammersley
if [ $? -ne 0 ]; then
  echo "Errors linking and loading hammersley_dataset.o"
  exit
fi
rm hammersley_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/hammersley_dataset
#
echo "Executable installed as ~/bin/$ARCH/hammersley_dataset"
