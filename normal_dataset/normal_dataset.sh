#!/bin/bash
#
gfortran -c -g normal_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling normal_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran normal_dataset.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading normal_dataset.o"
  exit
fi
rm normal_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/normal_dataset
#
echo "Executable installed as ~/bin/$ARCH/normal_dataset"
