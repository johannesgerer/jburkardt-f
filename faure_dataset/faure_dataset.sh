#!/bin/bash
#
gfortran -c -g faure_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling faure_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran faure_dataset.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading faure_dataset.o"
  exit
fi
rm faure_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/faure_dataset
#
echo "Executable installed as ~/bin/$ARCH/faure_dataset"
