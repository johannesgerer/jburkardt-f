#!/bin/bash
#
gfortran -c -g lcvt_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling lcvt_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran lcvt_dataset.o
if [ $? -ne 0 ]; then
  echo "Errors linking and loading lcvt_dataset.o"
  exit
fi
rm lcvt_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/lcvt_dataset
#
echo "Executable installed as ~/bin/$aRCH/lcvt_dataset"
