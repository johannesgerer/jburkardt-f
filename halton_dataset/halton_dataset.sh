#!/bin/bash
#
gfortran -c -g halton_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling halton_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran halton_dataset.o -L$HOME/lib/$ARCH -lhalton
if [ $? -ne 0 ]; then
  echo "Errors linking and loading halton_dataset.o"
  exit
fi
rm halton_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/halton_dataset
#
echo "Program installed as ~/bin/$ARCH/halton_dataset"
