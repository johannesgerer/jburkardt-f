#!/bin/bash
#
gfortran -c -g uniform_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling uniform_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran uniform_dataset.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading uniform_dataset.o"
  exit
fi
rm uniform_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/uniform_dataset
#
echo "Executable installed as ~/bin/$ARCH/uniform_dataset"
