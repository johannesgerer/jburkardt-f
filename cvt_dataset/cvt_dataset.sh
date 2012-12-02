#!/bin/bash
#
gfortran -c -g cvt_dataset.f90 >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling cvt_dataset.f90"
  exit
fi
rm compiler.txt
#
gfortran cvt_dataset.o -L$HOME/lib/$ARCH -lcvt
if [ $? -ne 0 ]; then
  echo "Errors linking and loading cvt_dataset.o"
  exit
fi
rm cvt_dataset.o
#
chmod ugo+x a.out
mv a.out ~/bin/$ARCH/cvt_dataset
#
echo "Program installed as ~/bin/$ARCH/cvt_dataset"
